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
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{scale}{Scaling passthrough texture}
 * \order{3}
 * \parameters{
 *     \parameter{value}{\Spectrum\Or\Texture}{
 *       Specifies the spectrum or nested texture that should be scaled
 *     }
 *     \parameter{value}{\Float}{
 *       Specifies the scale value
 *     }
 * }
 *
 * This simple plugin wraps a nested texture plugin and multiplies its
 * contents by a user-specified value. This can be quite useful when a
 * texture is too dark or too bright. The plugin can also be used to adjust
 * the height of a bump map when using the \pluginref{bumpmap} plugin.
 *
 * \begin{xml}[caption=Scaling the contents of a bitmap texture]
 * <texture type="scale">
 *     <float name="scale" value="0.5"/>
 *
 *     <texture type="bitmap">
 *         <string name="filename" value="wood.jpg"/>
 *     </texture>
 * </texture>
 * \end{xml}
 */

class Blend : public Texture {
public:

	Blend(const Properties &props) : Texture(props) {
		m_color0 = new ConstantSpectrumTexture(props.getSpectrum("color0", Spectrum(.5f)));
		m_color1 = new ConstantSpectrumTexture(props.getSpectrum("color1", Spectrum(.5f)));
		m_amount = props.getFloat("amount", 0.5);
	}

	Blend(Stream *stream, InstanceManager *manager)
		: Texture(stream, manager) {
		m_color0 = static_cast<Texture *>(manager->getInstance(stream));
		m_color1 = static_cast<Texture *>(manager->getInstance(stream));
		m_amount = stream->readFloat();
	}

	void configure() {
		if (m_color0 == NULL || m_color1 == NULL)
			Log(EError, "Blend plugin needs 2 textures!");
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if(name=="color0")
				m_color0 = static_cast<Texture *>(child);
			else if(name=="color1")
				m_color1 = static_cast<Texture *>(child);
		}
		else
			Texture::addChild(name, child);
	}
	
	Spectrum calculate(const Spectrum& color0, const Spectrum& color1) const{
		return color0 * (1.0f - m_amount) + color1 * m_amount;
	}

	Spectrum eval(const Intersection &its, bool filter) const {
		Spectrum color0 = m_color0->eval(its, filter);
		Spectrum color1 = m_color1->eval(its, filter);
		return  calculate(color0,color1);
	}

	void evalGradient(const Intersection &its, Spectrum *gradient) const {
		// m_nested->evalGradient(its, gradient);
		// gradient[0] *= m_scale;
		// gradient[1] *= m_scale;
	}

	Spectrum getAverage() const {
		Spectrum color0 = m_color0->getAverage();
		Spectrum color1 = m_color1->getAverage();
		return calculate(color0,color1);
	}

	Spectrum getMaximum() const {
		Spectrum color0 = m_color0->getMaximum();
		Spectrum color1 = m_color1->getMaximum();
		return calculate(color0,color1);
	}

	Spectrum getMinimum() const {
		Spectrum color0 = m_color0->getMinimum();
		Spectrum color1 = m_color1->getMinimum();
		return calculate(color0,color1);
	}

	bool isConstant() const {
		return m_color0->isConstant() && m_color1->isConstant();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Blend[" << endl
			<< "  color0 = " << indent(m_color0->toString()) << "," << endl
			<< "  color1 = " << indent(m_color1->toString()) << "," << endl
			<< "  amount = " << m_amount << endl
			<< "]";
		return oss.str();
	}

	bool usesRayDifferentials() const {
		return m_color0->usesRayDifferentials() || m_color1->usesRayDifferentials();
	}

	bool isMonochromatic() const {
		return m_color0->isMonochromatic() &&
			m_color1->isMonochromatic();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		manager->serialize(stream, m_color0.get());
		manager->serialize(stream, m_color1.get());
		stream->writeFloat(m_amount);
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ref<const Texture> m_color0, m_color1;
	Float m_amount;
};

// ================ Hardware shader implementation ================

class BlendShader : public Shader {
public:
	BlendShader(Renderer *renderer, const Texture *color0,
		const Texture *color1, Float amount) : Shader(renderer, ETextureShader),
		m_color0(color0), m_color1(color1), m_amount(amount){
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
			<< "uniform float " << evalName << "_amount;" << endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return " << depNames[0] << "(uv) * (1 - "<<evalName << "_amount) + "<< endl
			<< depNames[1] << "(uv) * "<<evalName << "_amount;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_amount", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
		int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_amount);
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_color0;
	ref<const Texture> m_color1;
	ref<Shader> m_color0Shader;
	ref<Shader> m_color1Shader;
	Float m_amount;
};

Shader *Blend::createShader(Renderer *renderer) const {
	return new BlendShader(renderer, m_color0.get(), m_color1.get(),
		m_amount);
}

MTS_IMPLEMENT_CLASS(BlendShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Blend, false, Texture2D)
MTS_EXPORT_PLUGIN(Blend, "Blend texture");
MTS_NAMESPACE_END
