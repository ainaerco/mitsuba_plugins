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
#include <mitsuba/core/spline.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/basicshader.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

MTS_NAMESPACE_BEGIN

/*!\plugin{ColorRamp}{ColorRamp}
 * \order{2}
 * \parameters{
 *     \parameter{color0, color1, color2, ...}{\Spectrum\Or\Texture}{
 *       Color values for ramp, should be named in proper order 0,1,2,...
 *     }
 *     \parameter{position0, position1, position2, ...}{\Float\Or\Texture}{
 *       Positions of colors in the ramp, should be named in proper order 0,1,2,...
 *     }
 *     \parameter{direction}{\String}{
 *       Direction or ramp
 *       \begin{enumerate}[(i)]
 *           \item \code{u}: Orient ramp in u direction
 *           \vspace{-1.5mm}
 *           \item \code{v}: Orient ramp in v direction
 *           \vspace{-1.5mm}
 *       \end{enumerate} u\default{0.5}
 *     }
 *     \parameter{input}{\Float\Or\Texture}{
 *       Position of lookup point
 *     }
 * }
 * \renderings{
 *     \rendering{ColorRamp applied to the material pf test object}{tex_ramp}
 * }
 * This plugin implements a simple ColorRamp texture with
 * connectable textures to colors, positions and input.
 */
class ColorRamp : public Texture {
public:
	ColorRamp(const Properties &props) : Texture(props) {
		m_numSamples = 0;
		std::string colorProp = (boost::format("color%d") % m_numSamples).str();
		std::string posProp = (boost::format("position%d") % m_numSamples).str();
		while (props.hasProperty(colorProp) && props.hasProperty(posProp))
		{
			m_numSamples++;
			m_colors.push_back(new ConstantSpectrumTexture(props.getSpectrum(colorProp)));
			m_positions.push_back(new ConstantFloatTexture(props.getFloat(posProp)));
			colorProp = (boost::format("color%d") % m_numSamples).str();
			posProp = (boost::format("position%d") % m_numSamples).str();
		}
		m_input = new ConstantFloatTexture(props.getFloat("input", 0.5f));
	}

	ColorRamp(Stream *stream, InstanceManager *manager)
	 : Texture(stream, manager) {
	 	m_numSamples = (size_t)stream->readUInt();
		m_colors.resize(m_numSamples);
		m_positions.resize(m_numSamples);
		for (size_t i=0; i<m_numSamples; ++i) {
			m_colors[i] = static_cast<Texture *>(manager->getInstance(stream));
		}
		for (size_t i=0; i<m_numSamples; ++i) {
			m_positions[i] = static_cast<Texture *>(manager->getInstance(stream));
		}
	 	m_input = static_cast<Texture *>(manager->getInstance(stream));
	 	configure();
	}

	void configure()
	{
		if (m_numSamples < 3)
			Log(EError, "Color Ramp: expect at least two points!");
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		stream->writeUInt((uint32_t) m_numSamples);
		BOOST_FOREACH(ref<Texture> tex, m_colors)
			manager->serialize(stream, tex.get());
		BOOST_FOREACH(ref<Texture> tex, m_positions)
			manager->serialize(stream, tex.get());
		manager->serialize(stream, m_input.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			boost::regex colorExpr("color(\\d+)$");
			boost::regex posExpr("position(\\d+)$");
			boost::match_results<std::string::const_iterator> what;

			if(name == "input") {
				m_input = static_cast<Texture *>(child);
			} else if(boost::regex_match(name, what, colorExpr)) {
					int i = boost::lexical_cast<int>(what[1]);
					m_colors[i] = static_cast<Texture *>(child);
			} else if(boost::regex_match(name, what, posExpr)) {
					int i = boost::lexical_cast<int>(what[1]);
					m_positions[i] = static_cast<Texture *>(child);
			} else
				Texture::addChild(name, child);
		} else
			Texture::addChild(name, child);
	}
	
	Spectrum eval(const Intersection &its, bool filter) const {
		std::vector<Float> positions;
		std::vector<Spectrum> colors;
		BOOST_FOREACH(ref<Texture> t, m_positions)
		{
			positions.push_back(t->eval(its, filter).average());
		}
		BOOST_FOREACH(ref<Texture> t, m_colors)
		{
			colors.push_back(t->eval(its, filter));
		}

		Spectrum result;
		for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
			std::vector<Float> samples;
			BOOST_FOREACH(Spectrum& s, colors)
				samples.push_back(s[i]);
			result[i] = std::max(0.0f, std::min(1.0f, evalCubicInterp1DN(
				-its.uv.y, &positions[0], &samples[0], m_numSamples)));
		}
		return result;
	}

	bool usesRayDifferentials() const {
		bool usesRayDifferentials = false;
		BOOST_FOREACH(ref<Texture> t, m_colors)
			usesRayDifferentials |= t->usesRayDifferentials();
		BOOST_FOREACH(ref<Texture> t, m_positions)
			usesRayDifferentials |= t->usesRayDifferentials();
		return usesRayDifferentials;
	}

	Spectrum getMaximum() const {
		Spectrum max(0.0f);
		BOOST_FOREACH(ref<Texture> t, m_colors) {
			for (int i=0; i<SPECTRUM_SAMPLES; ++i)
				max[i] = std::max(max[i], t->getMaximum()[i]);
		}
		return max;
	}

	Spectrum getMinimum() const {
		Spectrum min(1.0f);
		BOOST_FOREACH(ref<Texture> t, m_colors) {
			for (int i=0; i<SPECTRUM_SAMPLES; ++i)
				min[i] = std::min(min[i], t->getMinimum()[i]);
		}
		return min;
	}

	Spectrum getAverage() const {
		Spectrum average(0.0f);
		BOOST_FOREACH(ref<Texture> t, m_colors)
			average += t->getAverage();
		average /= (float)m_numSamples;
		return average;
	}

	bool isConstant() const {
		return false;
	}

	bool isMonochromatic() const {
		bool isMonochromatic = true;
		BOOST_FOREACH(ref<Texture> t, m_colors)
			isMonochromatic &= t->isMonochromatic();
		return isMonochromatic;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "ColorRamp[" << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	std::vector<ref<Texture> > m_colors;
	std::vector<ref<Texture> > m_positions;
	ref<Texture> m_input;
	size_t m_numSamples;
};

// ================ Hardware shader implementation ================

class CheckerShader : public Shader {
public:
	CheckerShader(Renderer *renderer, const Texture *color0) : Shader(renderer, ETextureShader),
		m_color0(color0) {
		
		m_color0Shader = renderer->registerShaderForResource(m_color0.get());
	}

	bool isComplete() const {
		return m_color0Shader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_color0.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_color0Shader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_color0;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "        return " << depNames[0] << "(uv);" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
		int &textureUnitOffset) const {
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_color0;
	ref<Shader> m_color0Shader;
};

Shader *ColorRamp::createShader(Renderer *renderer) const {
	return new CheckerShader(renderer, m_colors[0].get());
}

MTS_IMPLEMENT_CLASS(CheckerShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(ColorRamp, false, Texture2D)
MTS_EXPORT_PLUGIN(ColorRamp, "ColorRamp texture");
MTS_NAMESPACE_END
