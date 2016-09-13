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

#include <mitsuba/render/texture.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{Checker}{Checker}
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
 *     \rendering{Checker applied to the material test object
 *                as well as the ground plane}{tex_Checker}
 * }
 * This plugin implements a simple procedural Checker texture with
 * customizable colors.
 */
class Checker : public Texture {
public:
	Checker(const Properties &props) : Texture(props) {
		m_color0 = new ConstantSpectrumTexture(props.getSpectrum("color0", Spectrum(.4f)));
		m_color1 = new ConstantSpectrumTexture(props.getSpectrum("color1", Spectrum(.2f)));
		m_uvOffset = Point2(
			props.getFloat("uoffset", 0.0f),
			props.getFloat("voffset", 0.0f)
		);
		Float uvscale = props.getFloat("uvscale", 1.0f);
		m_uvScale = Vector2(
			props.getFloat("uscale", uvscale),
			props.getFloat("vscale", uvscale)
		);
	}

	Checker(Stream *stream, InstanceManager *manager)
	 : Texture(stream, manager) {
	 	m_color0 = static_cast<Texture *>(manager->getInstance(stream));
	 	m_color1 = static_cast<Texture *>(manager->getInstance(stream));
	 	m_uvOffset = Point2(stream);
		m_uvScale = Vector2(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		manager->serialize(stream, m_color0.get());
		manager->serialize(stream, m_color1.get());
		m_uvOffset.serialize(stream);
		m_uvScale.serialize(stream);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "color0")
				m_color0 = static_cast<Texture *>(child);
			else if (name == "color1")
				m_color1 = static_cast<Texture *>(child);
			else
			Texture::addChild(name, child);
		} else
			Texture::addChild(name, child);
	}
	
	Spectrum eval(const Intersection &its, bool filter) const {
		Intersection itst = its;
		itst.uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;

		int x = 2*math::modulo((int) (itst.uv.x * 2), 2) - 1,
			y = 2*math::modulo((int) (itst.uv.y * 2), 2) - 1;

		if (x*y == 1)
			return m_color0->eval(its, filter);
		else
			return m_color1->eval(its, filter);
	}

	bool usesRayDifferentials() const {
		return m_color0->usesRayDifferentials() || m_color1->usesRayDifferentials();
	}

	Spectrum getMaximum() const {
		Spectrum max0 = m_color0->getMaximum();
		Spectrum max1 = m_color1->getMaximum();
		Spectrum max;
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			max[i] = std::max(max0[i], max1[i]);
		return max;
	}

	Spectrum getMinimum() const {
		Spectrum min0 = m_color0->getMinimum();
		Spectrum min1 = m_color1->getMinimum();
		Spectrum min;
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			min[i] = std::min(min0[i], min1[i]);
		return min;
	}

	Spectrum getAverage() const {
		Spectrum avg0 = m_color0->getAverage();
		Spectrum avg1 = m_color1->getAverage();
		return (avg0 + avg1) * 0.5f;
	}

	bool isConstant() const {
		return m_color0->isConstant() && m_color1->isConstant();
	}

	bool isMonochromatic() const {
		return m_color0->isMonochromatic() &&
			m_color1->isMonochromatic();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Checker[" << endl
			<< "    color1 = " << indent(m_color1.toString()) << "," << endl
			<< "    color0 = " << indent(m_color0.toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ref<Texture> m_color1;
	ref<Texture> m_color0;
	Point2 m_uvOffset;
	Vector2 m_uvScale;
};

// ================ Hardware shader implementation ================

class CheckerShader : public Shader {
public:
	CheckerShader(Renderer *renderer, const Texture *color0,
		const Texture *color1, const Point2 &uvOffset,
		const Vector2 &uvScale) : Shader(renderer, ETextureShader),
		m_color0(color0), m_color1(color1),
		m_uvOffset(uvOffset), m_uvScale(uvScale) {
		
		m_color0Shader = renderer->registerShaderForResource(m_color0.get());
		m_color1Shader = renderer->registerShaderForResource(m_color1.get());
	}

	bool isComplete() const {
		return m_color0Shader.get() != NULL &&
			   m_color1Shader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_color0.get());
		renderer->unregisterShaderForResource(m_color1.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_color0Shader.get());
		deps.push_back(m_color1Shader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_color0;" << endl
			<< "uniform vec3 " << evalName << "_color1;" << endl
			<< "uniform vec2 " << evalName << "_uvOffset;" << endl
			<< "uniform vec2 " << evalName << "_uvScale;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    vec2 uvs = vec2(" << endl
			<< "        uv.x * " << evalName << "_uvScale.x + " << evalName << "_uvOffset.x," << endl
			<< "        uv.y * " << evalName << "_uvScale.y + " << evalName << "_uvOffset.y);" << endl
			<< "    float x = 2*(mod(int(uvs.x*2), 2)) - 1, y = 2*(mod(int(uvs.y*2), 2)) - 1;" << endl
			<< "    if (x*y == 1)" << endl
			<< "        return " << depNames[0] << "(uv);" << endl
			<< "    else" << endl
			<< "        return " << depNames[1] << "(uv);" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_uvOffset", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_uvScale", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
		int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_uvOffset);
		program->setParameter(parameterIDs[1], m_uvScale);
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_color0;
	ref<const Texture> m_color1;
	ref<Shader> m_color0Shader;
	ref<Shader> m_color1Shader;
	Point2 m_uvOffset;
	Vector2 m_uvScale;
};

Shader *Checker::createShader(Renderer *renderer) const {
	return new CheckerShader(renderer, m_color0.get(), m_color1.get(),
		m_uvOffset, m_uvScale);
}

MTS_IMPLEMENT_CLASS(CheckerShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Checker, false, Texture2D)
MTS_EXPORT_PLUGIN(Checker, "Checker texture");
MTS_NAMESPACE_END
