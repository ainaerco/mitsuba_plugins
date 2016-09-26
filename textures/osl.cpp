/*
	This file is part of Mitsuba, a physically based rendering system.

	Copyright (c) 2007-2014 by Wenzel Jakob and others.

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
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/thread/tss.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/join.hpp>

#include <OSL/oslexec.h>
#include <OSL/oslcomp.h>
#include <OSL/oslquery.h>

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/basicshader.h>
#include "simplerend.h"

MTS_NAMESPACE_BEGIN

using namespace OSL;

static ShadingSystem *shadingsys = NULL;
static ErrorHandler errhandler;
static SimpleRenderer rend;
static OSL::Matrix44 Mshad;  // "shader" space to "common" space matrix
static OSL::Matrix44 Mobj;   // "object" space to "common" space matrix

class PerThreadContext
{
public:
	PerThreadContext()
	{
		thread_info = shadingsys->create_thread_info();
		ctx = shadingsys->get_context (thread_info);
	}
	~PerThreadContext()
	{
		shadingsys->release_context (ctx);
		shadingsys->destroy_thread_info(thread_info);
	}
	OSL::PerThreadInfo* thread_info;
	ShadingContext* ctx;
};
boost::thread_specific_ptr<PerThreadContext> g_threadctx;

/*!\plugin{OSLTexture}{OSLTexture}
 * \order{2}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Filename of osl or oso shader
 *       \default{color(u,v,0)}
 *     }
 *     \parameter{searchpath}{\String}{
 *       Sets shaders search path for OSL shading system
 *     }
 *     \parameter{output}{\String}{
 *       Sets output of shader
 *       \default{first output in osl}
 *     }
 * }
 * This plugin implements an OSL texture. Closures, tracing and transformations 
 * are not supported yet.
 */

inline void
setup_transformations (SimpleRenderer &rend, OSL::Matrix44 &Mshad,
					   OSL::Matrix44 &Mobj)
{
	Matrix44 M (1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
	rend.camera_params (M, ustring("perspective"), 90.0f,
						0.1f, 1000.0f, 64, 64);

}

inline void
test_group_attributes (ShaderGroup *group)
{
	int nt = 0;
	if (shadingsys->getattribute (group, "num_textures_needed", nt)) {
		std::cout << "Need " << nt << " textures:\n";
		ustring *tex = NULL;
		shadingsys->getattribute (group, "textures_needed",
								  TypeDesc::PTR, &tex);
		for (int i = 0; i < nt; ++i)
			std::cout << "    " << tex[i] << "\n";
		int unk = 0;
		shadingsys->getattribute (group, "unknown_textures_needed", unk);
		if (unk)
			std::cout << "    and unknown textures\n";
	}
	int nclosures = 0;
	if (shadingsys->getattribute (group, "num_closures_needed", nclosures)) {
		std::cout << "Need " << nclosures << " closures:\n";
		ustring *closures = NULL;
		shadingsys->getattribute (group, "closures_needed",
								  TypeDesc::PTR, &closures);
		for (int i = 0; i < nclosures; ++i)
			std::cout << "    " << closures[i] << "\n";
		int unk = 0;
		shadingsys->getattribute (group, "unknown_closures_needed", unk);
		if (unk)
			std::cout << "    and unknown closures\n";
	}
	int nglobals = 0;
	if (shadingsys->getattribute (group, "num_globals_needed", nglobals)) {
		std::cout << "Need " << nglobals << " globals:\n";
		ustring *globals = NULL;
		shadingsys->getattribute (group, "globals_needed",
								  TypeDesc::PTR, &globals);
		for (int i = 0; i < nglobals; ++i)
			std::cout << "    " << globals[i] << "\n";
	}
	int nuser = 0;
	if (shadingsys->getattribute (group, "num_userdata", nuser) && nuser) {
		std::cout << "Need " << nuser << " user data items:\n";
		ustring *userdata_names = NULL;
		TypeDesc *userdata_types = NULL;
		int *userdata_offsets = NULL;
		bool *userdata_derivs = NULL;
		shadingsys->getattribute (group, "userdata_names",
								  TypeDesc::PTR, &userdata_names);
		shadingsys->getattribute (group, "userdata_types",
								  TypeDesc::PTR, &userdata_types);
		shadingsys->getattribute (group, "userdata_offsets",
								  TypeDesc::PTR, &userdata_offsets);
		shadingsys->getattribute (group, "userdata_derivs",
								  TypeDesc::PTR, &userdata_derivs);
		DASSERT (userdata_names && userdata_types && userdata_offsets);
		for (int i = 0; i < nuser; ++i)
			std::cout << "    " << userdata_names[i] << ' '
					  << userdata_types[i] << "  offset="
					  << userdata_offsets[i] << " deriv="
					  << userdata_derivs[i] << "\n";
	}
	int nattr = 0;
	if (shadingsys->getattribute (group, "num_attributes_needed", nattr) && nattr) {
		std::cout << "Need " << nattr << " attributes:\n";
		ustring *names = NULL;
		ustring *scopes = NULL;
		shadingsys->getattribute (group, "attributes_needed",
								  TypeDesc::PTR, &names);
		shadingsys->getattribute (group, "attribute_scopes",
								  TypeDesc::PTR, &scopes);
		DASSERT (names && scopes);
		for (int i = 0; i < nattr; ++i)
			std::cout << "    " << names[i] << ' '
					  << scopes[i] << "\n";

		int unk = 0;
		shadingsys->getattribute (group, "unknown_attributes_needed", unk);
		if (unk)
			std::cout << "    and unknown attributes\n";
	}
	int raytype_queries = 0;
	shadingsys->getattribute (group, "raytype_queries", raytype_queries);
	std::cout << "raytype() query mask: " << raytype_queries << "\n";
}

class OSLTexture : public Texture {
public:
	OSLTexture(const Properties &props) : Texture(props) {
		shadingsys = new ShadingSystem (&rend, NULL, &errhandler);
		register_closures(shadingsys);
		shadingsys->attribute("lockgeom", 1);
		if(props.hasProperty("searchpath")){
			shadingsys->attribute("searchpath:shader", 
				props.getString("searchpath"));
		}
		
		bool doCompile;
		std::string shadername;
		std::string sourcecode;
		m_outputVar = props.getString("output","result");
		fs::path filepath = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		if(fs::exists(filepath) && filepath.extension()==".osl") {
			shadername = filepath.stem().string();
			std::ifstream in (filepath.string().c_str(), 
							  std::ios::in | std::ios::binary);
			if (in) {
				std::ostringstream contents;
				contents << in.rdbuf();
				in.close ();
				sourcecode = contents.str();
			}
			doCompile = true;
		}
		else if(fs::exists(filepath) && filepath.extension()==".oso") {
			shadername = filepath.stem().string();
			doCompile = false;
		}
		else
		{
			shadername = "expr";
			sourcecode =
			"shader " + shadername + " (\n"
			"    float s = u [[ int lockgeom=0 ]],\n"
			"    float t = v [[ int lockgeom=0 ]],\n"
			"    output color " + m_outputVar + " = 0,\n"
			"  )\n"
			"{\n"
			"    " + m_outputVar + " = color(u,v,0);\n"
			"}\n";
			doCompile = true;
		}

		if(doCompile) {
			std::string osobuffer;
			OSLCompiler compiler;
			std::vector<std::string> options;
			std::string OSLHOME = getenv ("OSLHOME");
			if(!OSLHOME.empty())
				options.push_back("-I"+OSLHOME+"/shaders");
			if(!compiler.compile_buffer (sourcecode, osobuffer, options)) {
				Log(EError, "Could not compile %s", shadername);
			}
			if(!shadingsys->LoadMemoryCompiledShader (shadername, osobuffer)) {
				Log(EError, "Could not load compiled buffer from %s", shadername);
			}
		}
		m_shadergroup = shadingsys->ShaderGroupBegin();
		shadingsys->Shader ("surface", shadername);
		shadingsys->ShaderGroupEnd();
		// setup_transformations(rend, Mshad, Mobj);
	}

	OSLTexture(Stream *stream, InstanceManager *manager)
	 : Texture(stream, manager) {
		Log(EInfo, "OSL Stream Constructor not implemented!");
	}

	~OSLTexture()
	{
		m_shadergroup.reset();
		delete shadingsys;
	}

	void configure()
	{
		// Check output
		bool output_found = false;
		std::vector<std::string> outputs;
		OSLQuery q;
		q.init (m_shadergroup.get(), 0); // layer 0
		for (size_t p = 0;  p < q.nparams(); ++p) {
			const OSLQuery::Parameter *param = q.getparam(p);
			if(!param->isoutput)
				continue;
			if(m_outputVar == param->name) {
				output_found = true;
				break;
			}
			outputs.push_back(param->name.string());
		}
		if(!output_found) {
			if(outputs.size() == 0)
				Log(EError, "Outputs not found in OSL!"); 
			else
				m_outputVar = outputs[0];
			else
				Log(EError, "Wrong 'output' specified, available options are: %s", boost::algorithm::join(outputs, ", "));
		}

		std::vector<const char *> aovnames;
		aovnames.push_back(m_outputVar.c_str());
		shadingsys->attribute (NULL,"renderer_outputs",
							   TypeDesc(TypeDesc::STRING,(int)aovnames.size()),
							   &aovnames[0]);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		Texture::addChild(name, child);
	}
	
	Spectrum eval(const Intersection &its, bool filter) const {
		Spectrum result(0.0f);
		PerThreadContext* threadCtx = g_threadctx.get();
		if( !threadCtx ) {
			threadCtx = new PerThreadContext();
			g_threadctx.reset(threadCtx);
		}
		ShadingContext *ctx = threadCtx->ctx;

		ShaderGlobals sg;
		memset(&sg, 0, sizeof(ShaderGlobals));
		sg.P = Vec3(its.p.x, its.p.y, its.p.z);
		// sg.dPdx = Vec3(its.dpdu.x, its.dpdu.y, its.dpdu.z);
		// sg.dPdy = Vec3(its.dpdv.x, its.dpdv.y, its.dpdv.z);
		sg.u = its.uv.x;
		sg.v = 1.0f - its.uv.y;
		// sg.N = Vec3(its.shFrame.n.x, its.shFrame.n.y, its.shFrame.n.z);
		// sg.Ng = Vec3(its.geoFrame.n.x, its.geoFrame.n.y, its.geoFrame.n.z);
		sg.dudx = its.dudx;
		sg.dudy = its.dudy;
		sg.dvdx = its.dvdx;
		sg.dvdy = its.dvdy;
		// sg.I = Vec3(its.wi.x, its.wi.y, its.wi.z);
		// sg.dIdx = Vec3(its.shFrame.s.x, its.shFrame.s.y, its.shFrame.s.z);
		// sg.dIdy = Vec3(its.shFrame.t.x, its.shFrame.t.y, its.shFrame.t.z);
		// if (sg.backfacing) {
		// 	sg.N = -sg.N;
		// 	sg.Ng = -sg.Ng;
		// }
		// sg.raytype = 0;
		sg.renderstate = &sg;
		// const Triangle& tri = its.shape->getTriangles()[its.primIndex];
		// sg.shader2common = OSL::TransformationPtr (&Mshad);
		// sg.object2common = OSL::TransformationPtr (&Mobj);
		// sg.surfacearea = 1.0f; // use heron's formula  
		// sg.flipHandedness = false;
		shadingsys->execute (ctx, *m_shadergroup, sg);
		TypeDesc t;
		const void *data = shadingsys->get_symbol (*ctx, ustring(m_outputVar.c_str()), t);
			// test_group_attributes (m_shadergroup.get());
			// Log(EError, "Output result not found.");
		if (!data) {
			Log(EError, "Output result not found.");
		} else
		{
			const float* fdata = (const float *)data;
			result.fromLinearRGB(fdata[0], fdata[1], fdata[2]);
		}
		return result;
	}

	bool usesRayDifferentials() const {
		return true;
	}

	Spectrum getMaximum() const {
		return Spectrum(1.0f);
	}

	Spectrum getMinimum() const {
		return Spectrum(0.0f);
	}

	Spectrum getAverage() const {
		return Spectrum(0.5f);
	}

	bool isConstant() const {
		return false;
	}

	bool isMonochromatic() const {
		return false;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "OSL[" << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ShaderGroupRef m_shadergroup;
	std::string m_outputVar;
};

// ================ Hardware shader implementation ================

class OSLTextureShader : public Shader {
public:
	OSLTextureShader(Renderer *renderer) : Shader(renderer, ETextureShader) {
	}

	bool isComplete() const {
		return true;
	}

	void cleanup(Renderer *renderer) {
	}

	void putDependencies(std::vector<Shader *> &deps) {
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "        return vec3(0.2);" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
		int &textureUnitOffset) const {
	}

	MTS_DECLARE_CLASS()
};

Shader *OSLTexture::createShader(Renderer *renderer) const {
	return new OSLTextureShader(renderer);
}

MTS_IMPLEMENT_CLASS(OSLTextureShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(OSLTexture, false, Texture2D)
MTS_EXPORT_PLUGIN(OSLTexture, "OSL texture");
MTS_NAMESPACE_END
