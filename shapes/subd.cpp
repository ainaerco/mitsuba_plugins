/*
	This file is part of Mitsuba, a physically based rendering system.

	Copyright (c) 2007-2012 by Wenzel Jakob and others.

	Mitsuba is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License Version 3
	as published by the Free Software Foundation.

	Mitsuba is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sensor.h>

#include <Alembic/Abc/All.h>
#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreFactory/All.h>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

// OSD includes
#define not !
#define and &&
#include <osd/mesh.h>
#include <far/primvarRefiner.h>
#include <osd/cpuVertexBuffer.h>
#include <far/topologyDescriptor.h>
#include <far/stencilTableFactory.h>

#include "threadsafeMap.h"


MTS_NAMESPACE_BEGIN

using namespace OpenSubdiv;
namespace Abc  = Alembic::Abc;
namespace AbcA = Alembic::AbcCoreAbstract;
namespace AbcF = Alembic::AbcCoreFactory;
namespace AbcG = Alembic::AbcGeom;

// Vertex container implementation.
struct SubdVertex {

	SubdVertex() { }

	SubdVertex(float xx, float yy, float zz) {
		x = xx;
		y = yy;
		z = zz;
	}

	SubdVertex(SubdVertex const & src) {
		x = src.x;
		y = src.y;
		z = src.z;
	}

	void Clear( void * =0 ) {
		x=y=z=0.0f;
	}

	void AddWithWeight(SubdVertex const & src, float weight) {
		x+=weight*src.x;
		y+=weight*src.y;
		z+=weight*src.z;
	}

	float x,y,z;
	SubdVertex operator-(const SubdVertex &v) const {
		return SubdVertex(x - v.x,
			y - v.y,
			z - v.z);
	}
};

inline SubdVertex cross(const SubdVertex &v1, const SubdVertex &v2) {
	return SubdVertex(
		(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)
	);
}

struct SubdVertexUV {

	SubdVertexUV() { }

	SubdVertexUV(SubdVertexUV const & src) {
		u = src.u;
		v = src.v;
	}

	void Clear( void * =0 ) {
		u=v=0.0f;
	}

	void AddWithWeight(SubdVertexUV const & src, float weight) {
		u+=weight*src.u;
		v+=weight*src.v;
	}

	float u, v;
};

struct Vertex {
	Point p;
	Normal n;
	Point2 uv;
};

/// For using vertices as keys in an associative structure
struct vertex_key_order : public
	std::binary_function<Vertex, Vertex, bool> {
public:
	bool operator()(const Vertex &v1, const Vertex &v2) const {
		if (v1.p.x < v2.p.x) return true;
		else if (v1.p.x > v2.p.x) return false;
		if (v1.p.y < v2.p.y) return true;
		else if (v1.p.y > v2.p.y) return false;
		if (v1.p.z < v2.p.z) return true;
		else if (v1.p.z > v2.p.z) return false;
		if (v1.n.x < v2.n.x) return true;
		else if (v1.n.x > v2.n.x) return false;
		if (v1.n.y < v2.n.y) return true;
		else if (v1.n.y > v2.n.y) return false;
		if (v1.n.z < v2.n.z) return true;
		else if (v1.n.z > v2.n.z) return false;
		if (v1.uv.x < v2.uv.x) return true;
		else if (v1.uv.x > v2.uv.x) return false;
		if (v1.uv.y < v2.uv.y) return true;
		else if (v1.uv.y > v2.uv.y) return false;
		return false;
	}
};

struct OBJTriangle {
	int p[3];
	int n[3];
	int uv[3];

	inline OBJTriangle() {
		memset(this, 0, sizeof(OBJTriangle));
	}
};

namespace FVarChannels
{
	enum Channels
	{
		k_uv,
		k_normal,
		k_color,
		k_nchannels
	};
}

template <typename T> T convertMatrix(const Abc::M44d& m)
{
	T r;
	for(int i = 0; i < 4; ++i)
	{
		for(int j = 0; j < 4; ++j)
		{
			r.m[i][j] = (Float)m[j][i];
		}
	}
	return r;
}


struct AlembicPolymesh
{
	std::string identifier;
	AbcG::IPolyMeshSchema schema;
	Transform transform;
};

struct WorkerData
{
	Abc::P3fArraySamplePtr samplePos;
	Abc::Int32ArraySamplePtr sampleCounts;
	Abc::Int32ArraySamplePtr sampleIndices;

	bool hasTexcoords;
	Abc::V2fArraySamplePtr sampleUV; 
	Abc::UInt32ArraySamplePtr sampleUVIndices;

	bool hasNormals;
	Abc::N3fArraySamplePtr sampleNormal;
	Abc::UInt32ArraySamplePtr sampleNormalIndices;

	bool hasCreases;
	Abc::FloatArraySamplePtr creaseSharpnesses;
	Abc::Int32ArraySamplePtr creaseIndices;
	Abc::Int32ArraySamplePtr creaseLengths;

	TriMesh* mesh;
	std::string shader;
};

class Subdiv : public Shape {
public:
	void recurseObjectChildren(const Abc::IObject& obj, 
		AbcA::index_t sampleIndex,
		Transform t, std::vector<AlembicPolymesh>& map)
	{
		for (size_t i = 0; i < obj.getNumChildren(); ++i)
		{
			Transform childT = t;
			Abc::IObject child = obj.getChild(i);
			const std::string fullName = child.getFullName();
			if(AbcG::IXform::matches(child.getMetaData()))
			{
				AbcG::IXform xformObj = AbcG::IXform(child, Abc::kWrapExisting);
				AbcG::IXformSchema schema = xformObj.getSchema();
				AbcG::XformSample sample;
				schema.get(sample, sampleIndex);
				childT = childT * convertMatrix<Matrix4x4>(sample.getMatrix());
			}
			else if(AbcG::IPolyMesh::matches(child.getMetaData()))
			{
				AbcG::IPolyMesh polyObj(child, Abc::kWrapExisting);
				AlembicPolymesh polymesh;
				polymesh.schema = polyObj.getSchema();
				polymesh.identifier = child.getFullName();
				polymesh.transform = t;
				map.push_back(polymesh);
				continue;
			}
			recurseObjectChildren(child, sampleIndex, childT, map);
		}
	}

	void findObject(Abc::IObject object, const std::string& identifier,
		AbcA::index_t sampleIndex,
		Transform t, std::vector<AlembicPolymesh>& map)
	{
		std::vector<std::string> parts;
		boost::split(parts, identifier, boost::is_any_of("/\\"));
		for(size_t i = 1; i < parts.size(); ++i)
		{
			Alembic::Abc::IObject child(object, parts[i]);
			object = child;
			
			if(!object)
				//Log(EError, "Can't find identifier %s!", identifier);
				return;
			if(AbcG::IXform::matches(object.getMetaData()))
			{
				AbcG::IXform xformObj = AbcG::IXform(object, Abc::kWrapExisting);
				AbcG::IXformSchema schema = xformObj.getSchema();
				AbcG::XformSample sample;
				schema.get(sample, sampleIndex);
				t = t * convertMatrix<Matrix4x4>(sample.getMatrix());

			}
			else if(AbcG::IPolyMesh::matches(object.getMetaData()))
			{
				AbcG::IPolyMesh polyObj(object, Abc::kWrapExisting);
				AlembicPolymesh polymesh;
				polymesh.schema = polyObj.getSchema();
				polymesh.identifier = object.getFullName();
				polymesh.transform = t;
				map.push_back(polymesh);
				return;
			}
		}
	}

	void workerFuncSimple(
		int threadId,
		const Transform& objectToWorld,
		const std::string& identifier,
		WorkerData* workerData
		
		// const WorkerData& workerData
		)
	{
		int nfaces   = (int)workerData->sampleCounts->size();
		int nverts   = (int)workerData->samplePos->size();
		int nuvs     = (int)workerData->sampleUV->size();
		int nnormals = (int)workerData->sampleNormal->size();

		// Log(EInfo, "\tLoaded mesh: nfaces %d, nverts %d, uvs %d, normals %d", nfaces, nverts, nuvs, nnormals);

		std::vector<int> normalIndices;
		if(workerData->hasNormals)
		{
			if(workerData->sampleNormalIndices==NULL) {
				int id = 0;
				for( int i = 0; i < workerData->sampleCounts->size(); ++i ) 
				{
					for( int j = 0; j < (*workerData->sampleCounts)[i]; ++j ) 
					{
						normalIndices.push_back( id );
						id++;
					}
				}
			}
			else {
				for( int i = 0; i < workerData->sampleNormalIndices->size(); ++i ) 
				{
					normalIndices.push_back(workerData->sampleNormalIndices->get()[i]);
				}
			}
		}

		std::vector<Vertex> vertexBuffer;
		typedef std::map<Vertex, uint32_t, vertex_key_order> VertexMapType;
		VertexMapType vertexMap;
		AABB aabb;
		int numTriFaces = 0;
		int currentIndex = 0;
		std::vector<int> traingleIndices;

		for (int face = 0; face < nfaces; ++face) {
			int vertsPerFace = workerData->sampleCounts->get()[face];
			if(vertsPerFace<=4) {
				currentIndex += vertsPerFace;
				continue;
			}
			// std::cout<<vertsPerFace<<"\n";
			for(int i=0;i<vertsPerFace-2;++i) {
				uint32_t key;
				// std::cout<<"\t"<<i<<"\n";
				for(int j = 0; j<3; ++j) {
					int id;
					if(i==0)
						id = j;
					else if(i == vertsPerFace-3 && j == 2)
						id = 0;
					else
						id = j+i*2;
					// std::cout<<"\t\t"<<id<<"    "<<workerData->sampleIndices->get()[currentIndex+id]<<"\n";

					Vertex vertex;
					const Imath::Vec3<float>& v = workerData->samplePos->get()
						[workerData->sampleIndices->get()[currentIndex+id]];
					vertex.p.x = v.x;
					vertex.p.y = v.y;
					vertex.p.z = v.z;
					vertex.p = objectToWorld(vertex.p);
					aabb.expandBy(vertex.p);
					if(workerData->hasTexcoords) {
						const Imath::Vec2<float>& uv = workerData->sampleUV->get()
							[workerData->sampleUVIndices->get()[currentIndex+id]];
						vertex.uv.x = uv.x;
						vertex.uv.y = -uv.y;
					}
					if(workerData->hasNormals) {
						const Imath::Vec3<float>& n = workerData->sampleNormal->get()
							[normalIndices[currentIndex+id]];
						vertex.n.x = n.x;
						vertex.n.y = n.y;
						vertex.n.z = n.z;
						vertex.n = objectToWorld(vertex.n);
					}

					VertexMapType::iterator it = vertexMap.find(vertex);
					if (it != vertexMap.end()) {
						key = it->second;
					} else {
						key = (uint32_t) vertexBuffer.size();
						vertexMap[vertex] = key;
						vertexBuffer.push_back(vertex);
					}
					traingleIndices.push_back(key);
				}
				numTriFaces++;
			}
			currentIndex += vertsPerFace;
			
		}

		ref<TriMesh> mesh = new TriMesh(identifier,
			numTriFaces, numTriFaces * 3,
			workerData->hasNormals, workerData->hasTexcoords, false,
			!workerData->hasNormals, false);

		Point    *target_positions = mesh->getVertexPositions();
		Normal   *target_normals   = mesh->getVertexNormals();
		Point2   *target_texcoords = mesh->getVertexTexcoords();
		Triangle *target_triangles = mesh->getTriangles();

		mesh->getAABB() = aabb;
		for(int i=0;i<numTriFaces;i++) {
			Triangle& t = target_triangles[i];
			t.idx[0] = traingleIndices[i*3];
			t.idx[1] = traingleIndices[i*3+1];
			t.idx[2] = traingleIndices[i*3+2];
		}
		BOOST_FOREACH(Vertex& v,vertexBuffer) {
			*target_positions++ = v.p;
		}
		if(workerData->hasNormals)
			BOOST_FOREACH(Vertex& v,vertexBuffer) {
				*target_normals++ = v.n;
			}
		if(workerData->hasTexcoords)
			BOOST_FOREACH(Vertex& v,vertexBuffer) {
				*target_texcoords++ = v.uv;
			}

		mesh->incRef();
		workerData->mesh = mesh;
		m_cache.insert(threadId, workerData);
		m_numTriangles += numTriFaces;
		m_numVertices += numTriFaces * 3;
		m_numVerticesLow += nverts;
	}

	void workerFunc(
		int threadId,
		const int maxlevel,
		const Transform& objectToWorld,
		const std::string& identifier,
		WorkerData* workerData
		
		// const WorkerData& workerData
		)
	{
		// OpenSubdiv
		Sdc::SchemeType type = Sdc::SCHEME_CATMARK;

		Sdc::Options options;
		options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);
		options.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_NONE); // Sdc::FVAR_LINEAR_ALL
		options.SetCreasingMethod(Sdc::Options::CREASE_CHAIKIN);
		// options.SetTriangleSubdivision(Sdc::Options::TRI_SUB_SMOOTH);
		typedef Far::TopologyDescriptor Descriptor;
		Descriptor desc;
		desc.numVertices        = (int)workerData->samplePos->size();
		desc.numFaces           = (int)workerData->sampleCounts->size();
		desc.numVertsPerFace    = workerData->sampleCounts->get();
		desc.vertIndicesPerFace = workerData->sampleIndices->get();
		int numFVarChannels = 0;
		int uvChannelId     = -1;
		int normalChannelId = -1;
		if(workerData->hasCreases)
		{
			desc.numCreases = (int)workerData->creaseLengths->size();
			desc.creaseVertexIndexPairs = workerData->creaseIndices->get();
			desc.creaseWeights = workerData->creaseSharpnesses->get();
		}

		std::vector<Descriptor::FVarChannel> channels;
		
		// UVs FVar
		if(workerData->hasTexcoords)
		{
			Descriptor::FVarChannel channel;
			channel.valueIndices = (Far::Index *)workerData->sampleUVIndices->get();
			channel.numValues = (int)workerData->sampleUV->size();
			channels.push_back(channel);
			uvChannelId = numFVarChannels;
			numFVarChannels++;
		}

		// Normals FVar
		std::vector<int> normalIndices;
		if(workerData->hasNormals)
		{
			Descriptor::FVarChannel channel;
			if(workerData->sampleNormalIndices==NULL) {
				int id = 0;
				for( int i = 0; i < workerData->sampleCounts->size(); ++i ) 
				{
					for( int j = 0; j < (*workerData->sampleCounts)[i]; ++j ) 
					{
						normalIndices.push_back( id );
						id++;
					}
				}
				channel.valueIndices = &normalIndices[0];
			}
			else {
				channel.valueIndices = (Far::Index *)workerData->sampleNormalIndices->get();
			}
			channel.numValues = (int)workerData->sampleNormal->size();
			channels.push_back(channel);
			normalChannelId = numFVarChannels;
			numFVarChannels++;
		}

		if(desc.numFVarChannels = numFVarChannels)
			desc.fvarChannels = &channels[0];

		Far::TopologyRefiner * refiner = 
			Far::TopologyRefinerFactory<Descriptor>::Create(desc,
				Far::TopologyRefinerFactory<Descriptor>::Options(type, options));

		// Uniformly refine the topolgy up to 'maxlevel'
		// note: fullTopologyInLastLevel must be true to work with face-varying data
		{
			Far::TopologyRefiner::UniformOptions refineOptions(maxlevel);
			refineOptions.fullTopologyInLastLevel = true;
			refiner->RefineUniform(refineOptions);
		}
		
		std::vector<SubdVertex> vbuffer(refiner->GetNumVerticesTotal());
		SubdVertex * verts = &vbuffer[0];
		std::memcpy(verts, workerData->samplePos->get(), sizeof(Abc::V3f) * desc.numVertices);

		// Allocate and initialize the UVs channel of 'face-varying' primvar data
		SubdVertexUV * fvVertsUV = NULL;
		std::vector<SubdVertexUV> fvBufferUV;
		if(workerData->hasTexcoords) {
			fvBufferUV.resize(refiner->GetNumFVarValuesTotal(uvChannelId));
			fvVertsUV = &fvBufferUV[0];
			std::memcpy(fvVertsUV, workerData->sampleUV->get(), sizeof(Abc::V2f) * desc.fvarChannels[uvChannelId].numValues);
		}

		// Allocate and initialize the N channel of 'face-varying' primvar data
		SubdVertex * fvVertsN = NULL;
		std::vector<SubdVertex> fvBufferN;
		if(workerData->hasNormals) {
			fvBufferN.resize(refiner->GetNumFVarValuesTotal(normalChannelId));
			fvVertsN = &fvBufferN[0];
			std::memcpy(fvVertsN, workerData->sampleNormal->get(), sizeof(Abc::V3f) * desc.fvarChannels[normalChannelId].numValues);
		}
		// Interpolate vertex and face-varying primvar data
		Far::PrimvarRefiner primvarRefiner(*refiner);

		SubdVertex* srcVert = verts;
		SubdVertexUV* srcFVarUV = fvVertsUV;
		SubdVertex* srcFVarN = fvVertsN;
		for (int level = 1; level <= maxlevel; ++level) 
		{
			SubdVertex*     dstVert = srcVert + refiner->GetLevel(level-1).GetNumVertices();
			primvarRefiner.Interpolate(level, srcVert, dstVert);
			srcVert   = dstVert;
			
			if(srcFVarUV != NULL) {
				SubdVertexUV* dstFVarUV = srcFVarUV + refiner->GetLevel(level-1).GetNumFVarValues(uvChannelId);
				primvarRefiner.InterpolateFaceVarying(level, srcFVarUV, dstFVarUV, uvChannelId);
				srcFVarUV = dstFVarUV;
			}
			if(srcFVarN != NULL) {
				SubdVertex* dstFVarN    = srcFVarN + refiner->GetLevel(level-1).GetNumFVarValues(normalChannelId);
				primvarRefiner.InterpolateFaceVarying(level, srcFVarN, dstFVarN, normalChannelId);
				srcFVarN  = dstFVarN;
			}
		}

		Far::TopologyLevel const & refLastLevel = refiner->GetLevel(maxlevel);
		int nfaces, nverts, nuvs, nnormals, firstOfLastVert, firstOfLastUvs, firstOfLastN;

		nfaces              = refLastLevel.GetNumFaces();
		nverts              = refLastLevel.GetNumVertices();
		firstOfLastVert     = refiner->GetNumVerticesTotal() - nverts;
		if(workerData->hasTexcoords) {
			nuvs            = refLastLevel.GetNumFVarValues(uvChannelId);
			firstOfLastUvs  = refiner->GetNumFVarValuesTotal(uvChannelId) - nuvs;
		}

		if(workerData->hasNormals) {
			nnormals        = refLastLevel.GetNumFVarValues(normalChannelId);
			firstOfLastN    = refiner->GetNumFVarValuesTotal(normalChannelId) - nnormals;
		}

		// Limit position, derivatives and normals
	    std::vector<SubdVertex> fineLimitPos(nverts);
	    std::vector<SubdVertex> fineDu(nverts);
	    std::vector<SubdVertex> fineDv(nverts);
	    std::vector<SubdVertex> finePosBuffer(nverts);
	    primvarRefiner.Interpolate( maxlevel, srcVert, finePosBuffer);
	    primvarRefiner.Limit(finePosBuffer, fineLimitPos, fineDu, fineDv);

		// Log(EInfo, "\tLoaded mesh: nfaces %d, nverts %d, uvs %d, normals %d", nfaces, nverts, nuvs, nnormals);
		ref<TriMesh> mesh = new TriMesh(identifier,
			nfaces * 2, nfaces * 2 * 3,
			workerData->hasNormals, workerData->hasTexcoords, false,
			!workerData->hasNormals, false);

		Point    *target_positions = mesh->getVertexPositions();
		Normal   *target_normals   = mesh->getVertexNormals();
		Point2   *target_texcoords = mesh->getVertexTexcoords();
		Triangle *target_triangles = mesh->getTriangles();

		AABB& aabb = mesh->getAABB();

		std::vector<Vertex> vertexBuffer;
		typedef std::map<Vertex, uint32_t, vertex_key_order> VertexMapType;
		VertexMapType vertexMap;

		for (int face = 0; face < nfaces; ++face) {
			Far::ConstIndexArray fverts   = refLastLevel.GetFaceVertices(face);
			Far::ConstIndexArray fuvs;
			Far::ConstIndexArray fnormals;
			if(workerData->hasTexcoords) 
				fuvs     = refLastLevel.GetFaceFVarValues(face, uvChannelId);
			if(workerData->hasNormals) 
				fnormals = refLastLevel.GetFaceFVarValues(face, normalChannelId);
			if(fverts.size()!=4 || fuvs.size()!=4)
				continue;
			// std::cout<<"face "<<face<<"\n";
			for(int i=0;i<2;++i){
			
				Triangle& t1 = target_triangles[face * 2 + i];
				uint32_t key;
				for(int j = 0; j<3; ++j)
				{
					int id = j*(1-i)+ i*((j+2)%4);
					Vertex vertex;
					const SubdVertex& v  = verts[firstOfLastVert + fverts[id]];
					vertex.p.x = v.x;
					vertex.p.y = v.y;
					vertex.p.z = v.z;
					vertex.p = objectToWorld(vertex.p);
					aabb.expandBy(vertex.p);
					if(workerData->hasTexcoords) {
						const SubdVertexUV& uv = fvVertsUV[firstOfLastUvs + fuvs[id]];
						vertex.uv.x = uv.u;
						vertex.uv.y = -uv.v;
					}
					if(workerData->hasNormals) {
						SubdVertex n  = cross(fineDu[fverts[id]],fineDv[fverts[id]]);

						// const SubdVertex& n  = fvVertsN[firstOfLastN + fnormals[id]];
						vertex.n.x = -n.x;
						vertex.n.y = -n.y;
						vertex.n.z = -n.z;
						vertex.n = objectToWorld(vertex.n);
					}

					VertexMapType::iterator it = vertexMap.find(vertex);
					if (it != vertexMap.end()) {
						key = it->second;
					} else {
						key = (uint32_t) vertexBuffer.size();
						vertexMap[vertex] = key;
						vertexBuffer.push_back(vertex);
					}

					t1.idx[j] = key;
				}
			}

			
		}
		BOOST_FOREACH(Vertex& v, vertexBuffer) {
			*target_positions++ = v.p;
		}
		if(workerData->hasNormals)
			BOOST_FOREACH(Vertex& v, vertexBuffer) {
				*target_normals++ = v.n;
			}
		if(workerData->hasTexcoords)
			BOOST_FOREACH(Vertex& v, vertexBuffer) {
				*target_texcoords++ = v.uv;
			}

		mesh->incRef();
		workerData->mesh = mesh;
		m_cache.insert(threadId, workerData);
		m_numTriangles += nfaces;
		m_numVertices += nfaces * 2 * 3;
		m_numVerticesLow += nverts;
		delete refiner;
	}

	Subdiv(const Properties &props) : Shape(props), m_numTriangles(0), m_numVertices(0), m_numVerticesLow(0) {
		
		fs::path filePath = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		if (!fs::exists(filePath))
			Log(EError, "Alembic file \"%s\" could not be found!", filePath.string().c_str());
		m_name = filePath.stem().string();

		// Object-space -> World-space transformation
		const Transform transform = props.getTransform("toWorld", Transform());

		const int maxlevel = props.getInteger("level", 0);

		const std::string identifiersStr = props.getString("identifiers","");

		ref<Timer> timer = new Timer();
		size_t trianglesTotal = 0;

		AbcA::index_t sampleIndex = 0;

		AbcF::IFactory factory;
		AbcF::IFactory::CoreType coreType;
		
		Abc::IArchive archive = factory.getArchive(filePath.string(), coreType);
		Alembic::Abc::IObject object = archive.getTop();

		std::vector<AlembicPolymesh> alembicPolymeshes;
		if(!identifiersStr.empty())
		{
			std::vector<std::string> identifiers;
			boost::split(identifiers, identifiersStr, boost::is_any_of(";"));
			BOOST_FOREACH( std::string& identifier, identifiers )
			{
				findObject(object, identifier, sampleIndex, transform, alembicPolymeshes);
			}
		}
		else
			recurseObjectChildren(object, sampleIndex, transform, alembicPolymeshes);
		
		int numthreads = 4;
		for(size_t i = 0; i < alembicPolymeshes.size() / numthreads + 1; ++i)
		{
			boost::thread_group m_worker_threads;
			m_cache.clear();

			// AbcG::IPolyMeshSchema polySchema[4];
			
			for (int j = 0; (i < alembicPolymeshes.size() / numthreads && j < numthreads) || 
				(i == alembicPolymeshes.size() / numthreads && j < alembicPolymeshes.size() % numthreads); ++j)
			{
				WorkerData* workerData = new WorkerData();
				AbcG::IPolyMeshSchema polySchema = alembicPolymeshes[i * numthreads + j].schema;
				workerData->samplePos = polySchema.getPositionsProperty().getValue(sampleIndex);
				workerData->sampleCounts = polySchema.getFaceCountsProperty().getValue(sampleIndex);
				workerData->sampleIndices = polySchema.getFaceIndicesProperty().getValue(sampleIndex);
				AbcG::IN3fGeomParam normalsParam = polySchema.getNormalsParam();
				AbcG::IV2fGeomParam uvsParam = polySchema.getUVsParam();
				Abc::ICompoundProperty arbProp = polySchema.getArbGeomParams();

				if(arbProp.valid())
				{
					std::size_t numProps = arbProp.getNumProperties();
					for (std::size_t propId = 0; propId < numProps; ++propId)
					{
						const Alembic::Abc::PropertyHeader & propHeader = arbProp.getPropertyHeader(propId);
						std::string propName = propHeader.getName();
						if(propName!="shader")
							continue;
						Abc::DataType type = propHeader.getDataType();
						Alembic::Abc::IArrayProperty iProp(arbProp, propName);
						Alembic::AbcCoreAbstract::ArraySamplePtr samp;
						iProp.get(samp);
						switch (type.getPod())
						{
							case Alembic::Util::kStringPOD:
							{
								Alembic::Util::string * strData = (Alembic::Util::string *) samp->getData();
								workerData->shader = std::string(*strData);
							} break;
							default: break;
						}
					}
				}

				workerData->hasNormals = normalsParam.valid() && normalsParam.getValueProperty().valid();
				workerData->hasTexcoords = uvsParam.valid() && uvsParam.getValueProperty().valid();
				workerData->hasCreases = polySchema.getPropertyHeader(".creaseIndices") != NULL &&
					polySchema.getPropertyHeader(".creaseLengths") != NULL &&
					polySchema.getPropertyHeader(".creaseSharpnesses") != NULL;

				if(workerData->hasCreases)
				{
					Abc::IInt32ArrayProperty creaseIndicesProperty(polySchema, ".creaseIndices");
					Abc::IInt32ArrayProperty creaseLengthsProperty(polySchema, ".creaseLengths");
					Abc::IFloatArrayProperty creaseSharpnessesProperty(polySchema, ".creaseSharpnesses");
					creaseSharpnessesProperty.get(workerData->creaseSharpnesses, sampleIndex);
					creaseLengthsProperty.get(workerData->creaseLengths, sampleIndex);
					creaseIndicesProperty.get(workerData->creaseIndices, sampleIndex);
				}

				if(workerData->hasTexcoords)
				{
					workerData->sampleUV = uvsParam.getValueProperty().getValue( sampleIndex );
					if(uvsParam.getIndexProperty().valid())
						workerData->sampleUVIndices = uvsParam.getIndexProperty().getValue( sampleIndex );
				}
				if(workerData->hasNormals)
				{
					workerData->sampleNormal = normalsParam.getValueProperty().getValue( sampleIndex );
					if(normalsParam.getIndexProperty().valid())
						workerData->sampleNormalIndices = normalsParam.getIndexProperty().getValue( sampleIndex );
				}
				if(maxlevel) {
					m_worker_threads.create_thread(boost::bind(&Subdiv::workerFunc, this, 
						j,
						maxlevel,
						alembicPolymeshes[i * numthreads + j].transform,
						alembicPolymeshes[i * numthreads + j].identifier,
						workerData
					));
				} else {
					m_worker_threads.create_thread(boost::bind(&Subdiv::workerFuncSimple, this, 
						j,
						alembicPolymeshes[i * numthreads + j].transform,
						alembicPolymeshes[i * numthreads + j].identifier,
						workerData
					));
				}
				// Subdiv::workerFuncSimple(
				// 		j,
				// 		alembicPolymeshes[i * numthreads + j].transform,
				// 		alembicPolymeshes[i * numthreads + j].identifier,
				// 		workerData);
			}

			m_worker_threads.join_all();
			
			for (int j = 0; j < numthreads; ++j)
			{
				WorkerData* w = m_cache.find(j);
				if(w == NULL || !w->mesh)
					continue;
				m_meshes.push_back(w->mesh);
				m_materialAssignment.push_back(w->shader);
				delete w;
			}

		}
		Log(EInfo, "Done with %d triangles (took %i ms)", m_numTriangles, timer->getMilliseconds());
		std::cout<<"Datasize triangles "<<m_numTriangles*sizeof(Triangle)/1000000<<
			" Mb  vertices "<<m_numVertices*sizeof(Point)/1000000<<" Mb  additional data "<<m_numVertices*sizeof(Point2)/1000000<<" Mb\n";
		std::cout<<"LowDatasize vertices "<<m_numVerticesLow*sizeof(Point)/1000000<<" Mb  additional data "<<m_numVerticesLow*sizeof(Point2)/1000000<<" Mb\n";
	}

	Subdiv(Stream *stream, InstanceManager *manager) : Shape(stream, manager) {
		m_aabb = AABB(stream);
		m_name = stream->readString();
		uint32_t meshCount = stream->readUInt();
		m_meshes.resize(meshCount);

		for (uint32_t i=0; i<meshCount; ++i) {
			m_meshes[i] = static_cast<TriMesh *>(manager->getInstance(stream));
			m_meshes[i]->incRef();
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);

		m_aabb.serialize(stream);
		stream->writeString(m_name);
		stream->writeUInt((uint32_t) m_meshes.size());
		for (size_t i=0; i<m_meshes.size(); ++i)
			manager->serialize(stream, m_meshes[i]);
	}

	virtual ~Subdiv() {
		for (size_t i=0; i<m_meshes.size(); ++i)
			m_meshes[i]->decRef();
	}

	void configure() {
		Shape::configure();

		m_aabb.reset();
		for (size_t i=0; i<m_meshes.size(); ++i) {
			m_meshes[i]->configure();
			m_aabb.expandBy(m_meshes[i]->getAABB());
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		addChild(name, child, true);
	}

	void addChild(const std::string &name, ConfigurableObject *child, bool warn) {
		/// TODO: review this function
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(MTS_CLASS(BSDF))) {
			Shape::addChild(name, child);
			Assert(m_meshes.size() > 0);
			if (name == "") {
				for (size_t i=0; i<m_meshes.size(); ++i)
					m_meshes[i]->addChild(name, child);
			} else {
				bool found = false;
				for (size_t i=0; i<m_meshes.size(); ++i) {
					if (m_materialAssignment[i] == name) {
						found = true;
						m_meshes[i]->addChild(name, child);
					}
				}
				if (!found && warn)
					Log(EWarn, "Attempted to register the material named "
						"'%s', which does not occur in the OBJ file!", name.c_str());
			}
			m_bsdf->setParent(NULL);
		} else if (cClass->derivesFrom(MTS_CLASS(Emitter))) {
			if (m_meshes.size() > 1)
				Log(EError, "Cannot attach an emitter to an OBJ file "
					"containing multiple objects!");
			m_emitter = static_cast<Emitter *>(child);
			child->setParent(m_meshes[0]);
			m_meshes[0]->addChild(name, child);
		} else if (cClass->derivesFrom(MTS_CLASS(Sensor))) {
			if (m_meshes.size() > 1)
				Log(EError, "Cannot attach an sensor to an OBJ file "
					"containing multiple objects!");
			m_sensor = static_cast<Sensor *>(child);
			child->setParent(m_meshes[0]);
			m_meshes[0]->addChild(name, child);
		} else if (cClass->derivesFrom(MTS_CLASS(Subsurface))) {
			Assert(m_subsurface == NULL);
			m_subsurface = static_cast<Subsurface *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) {
				child->setParent(m_meshes[i]);
				m_meshes[i]->addChild(name, child);
			}
		} else if (cClass->derivesFrom(MTS_CLASS(Medium))) {
			Shape::addChild(name, child);
			for (size_t i=0; i<m_meshes.size(); ++i)
				m_meshes[i]->addChild(name, child);
		} else {
			Shape::addChild(name, child);
		}
	}

	bool isCompound() const {
		return true;
	}

	Shape *getElement(int index) {
		if (index >= (int) m_meshes.size())
			return NULL;
		Shape *shape = m_meshes[index];
		BSDF *bsdf = shape->getBSDF();
		Emitter *emitter = shape->getEmitter();
		Subsurface *subsurface = shape->getSubsurface();
		if (bsdf)
			bsdf->setParent(shape);
		if (emitter)
			emitter->setParent(shape);
		if (subsurface)
			subsurface->setParent(shape);
		return shape;
	}

	std::string getName() const {
		return m_name;
	}

	AABB getAABB() const {
		return m_aabb;
	}

	Float getSurfaceArea() const {
		Float sa = 0;
		for (size_t i=0; i<m_meshes.size(); ++i)
			sa += m_meshes[i]->getSurfaceArea();
		return sa;
	}

	size_t getPrimitiveCount() const {
		size_t result = 0;
		for (size_t i=0; i<m_meshes.size(); ++i)
			result += m_meshes[i]->getPrimitiveCount();
		return result;
	}

	size_t getEffectivePrimitiveCount() const {
		size_t result = 0;
		for (size_t i=0; i<m_meshes.size(); ++i)
			result += m_meshes[i]->getEffectivePrimitiveCount();
		return result;
	}

	MTS_DECLARE_CLASS()
private:
	std::vector<TriMesh*> m_meshes;
	std::vector<std::string> m_materialAssignment;
	std::string m_name;
	AABB m_aabb;
	MapSafe<int, WorkerData*> m_cache;

	boost::atomic<size_t> m_numTriangles;
	boost::atomic<size_t> m_numVertices;
	boost::atomic<size_t> m_numVerticesLow;
};

MTS_IMPLEMENT_CLASS_S(Subdiv, false, Shape)
MTS_EXPORT_PLUGIN(Subdiv, "Subdiv mesh loader");
MTS_NAMESPACE_END
