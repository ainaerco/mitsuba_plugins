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
#include <boost/atomic.hpp>


MTS_NAMESPACE_BEGIN

namespace Abc  = Alembic::Abc;
namespace AbcA = Alembic::AbcCoreAbstract;
namespace AbcF = Alembic::AbcCoreFactory;
namespace AbcG = Alembic::AbcGeom;

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

class AlembicShape : public Shape {
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
		const Transform& objectToWorld,
		const std::string& identifier,
		WorkerData* workerData)
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
			for(int i=0;i<vertsPerFace-2;++i) {
				uint32_t key;
				for(int j = 0; j<3; ++j) {
					int id;
					if(i==0)
						id = j;
					else if(i == vertsPerFace-3 && j == 2)
						id = 0;
					else
						id = j+i*2;

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
		m_numTriangles += numTriFaces;
		m_numVertices += numTriFaces * 3;
		m_numVerticesLow += nverts;
	}

	AlembicShape(const Properties &props) : Shape(props), m_numTriangles(0), m_numVertices(0), m_numVerticesLow(0) {
		
		fs::path filePath = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		if (!fs::exists(filePath))
			Log(EError, "Alembic file \"%s\" could not be found!", filePath.string().c_str());
		m_name = filePath.stem().string();

		// Object-space -> World-space transformation
		const Transform transform = props.getTransform("toWorld", Transform());
		const std::string identifiersStr = props.getString("identifiers","");

		ref<Timer> timer = new Timer();

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

		for(size_t i = 0; i < alembicPolymeshes.size(); ++i)
		{
			WorkerData* workerData = new WorkerData();
			AbcG::IPolyMeshSchema polySchema = alembicPolymeshes[i].schema;
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
					if(propName != "shader")
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

			workerFuncSimple(alembicPolymeshes[i].transform, 
				alembicPolymeshes[i].identifier,
				workerData);
			m_meshes.push_back(workerData->mesh);
			m_materialAssignment.push_back(workerData->shader);
			
		}
		
		Log(EInfo, "Done with %d triangles (took %i ms)", m_numTriangles, timer->getMilliseconds());
		std::cout<<"Datasize triangles "<<m_numTriangles*sizeof(Triangle)/1000000<<
			" Mb  vertices "<<m_numVertices*sizeof(Point)/1000000<<" Mb  additional data "<<m_numVertices*sizeof(Point2)/1000000<<" Mb\n";
		std::cout<<"LowDatasize vertices "<<m_numVerticesLow*sizeof(Point)/1000000<<" Mb  additional data "<<m_numVerticesLow*sizeof(Point2)/1000000<<" Mb\n";
	}

	AlembicShape(Stream *stream, InstanceManager *manager) : Shape(stream, manager) {
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

	virtual ~AlembicShape() {
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

	boost::atomic<size_t> m_numTriangles;
	boost::atomic<size_t> m_numVertices;
	boost::atomic<size_t> m_numVerticesLow;
};

MTS_IMPLEMENT_CLASS_S(AlembicShape, false, Shape)
MTS_EXPORT_PLUGIN(AlembicShape, "AlembicShape mesh loader");
MTS_NAMESPACE_END
