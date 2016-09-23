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
#include <OSL/oslexec.h>
#include <OSL/oslexec_pvt.h>

#include <mitsuba/render/texture.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

using namespace OSL;

/*!\plugin{OSLTexture}{OSLTexture}
 * \order{2}
 * \parameters{
 *     \parameter{color0, color1}{\Spectrum}{
 *       Color values for the two differently-colored patches
 *       \default{0.4 and 0.2}
 *     }
 *     \parameter{uoffset, voffset}{\Float}{
 *       Numerical offset that should be applied to UV values before a lookup
 *     }
 *     \parameter{uscale, vscale}{\Float}{
 *       Multiplicative factors that should be applied to UV values before a lookup
 *     }
 * }
 * \renderings{
 *     \rendering{OSLTexture applied to the material test object
 *                as well as the ground plane}{tex_OSLTexture}
 * }
 * This plugin implements a simple procedural OSLTexture texture with
 * customizable colors.
 */

class OSLTexture : public Texture {
public:
	OSLTexture(const Properties &props) : Texture(props) {
		
	}

	OSLTexture(Stream *stream, InstanceManager *manager)
	 : Texture(stream, manager) {
	 	
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		Texture::addChild(name, child);
	}
	
	Spectrum eval(const Intersection &its, bool filter) const {
		ShaderGlobals sg_osl;
		memset(&sg_osl, 0, sizeof(ShaderGlobals));
		sg_osl.P = its.p;
		sg_osl.dPdx = its.dpdu;
		sg_osl.dPdy = its.dpdv;
		sg_osl.u = its.uv.x;
		sg_osl.v = its.uv.y;
		sg_osl.N = its.shFrame.n;
		sg_osl.Ng = its.geoFrame.n;
		sg_osl.dudx = its.dudx;
		sg_osl.dudy = its.dudy;
		sg_osl.dvdx = its.dvdx;
		sg_osl.I = its.wi;
		sg_osl.dIdx = its.shFrame.s;
		sg_osl.dIdy = its.shFrame.t;
		if (sg_osl.backfacing) {
			sg_osl.N = -sg_osl.N;
			sg_osl.Ng = -sg_osl.Ng;
		}
		sg_osl.renderstate = &sg_osl;
		const Triangle& tri = its.shape->getTriangles()[its.primIndex];
		//sg_osl.surfacearea  // use heron's formula  
		sg_osl.flipHandedness = false;
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
	std::String m_oslString;
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
			<< "        return " << depNames[1] << "(uv);" << endl
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
MTS_EXPORT_PLUGIN(OSLTexture, "OSLTexture texture");
MTS_NAMESPACE_END
