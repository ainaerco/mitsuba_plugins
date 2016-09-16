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
#include <mitsuba/core/random.h>
#include <mitsuba/hw/basicshader.h>
#include <boost/functional/hash.hpp>
#include <boost/thread/tss.hpp>
#include "hsv.h"

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
boost::thread_specific_ptr<std::map<std::string, float> > g_threadMap;

class VariationTexture : public Texture {
public:
    VariationTexture(const Properties &props) : Texture(props) {
        if (props.hasProperty("value"))
            m_nested = new ConstantSpectrumTexture(
                props.getSpectrum("value"));

        if (props.hasProperty("scale") && props.getType("scale") == Properties::EFloat)
            m_scale = Spectrum(props.getFloat("scale", 1.0f));
        else
            m_scale = props.getSpectrum("scale", Spectrum(1.0f));

        Assert(m_scale.min() > 0);
        m_gamma = props.getFloat("gamma", 0.0f);
        m_hue = props.getFloat("hue", 0.0f);
        m_saturation = props.getFloat("saturation", 0.0f);
        m_exposure = props.getFloat("exposure", 0.0f);
        m_offset = props.getFloat("offset", 0.0f);
    }

    VariationTexture(Stream *stream, InstanceManager *manager)
        : Texture(stream, manager) {
        m_nested = static_cast<Texture *>(manager->getInstance(stream));
        m_scale = Spectrum(stream);
        m_gamma = stream->readFloat();
        m_hue = stream->readFloat();
        m_saturation = stream->readFloat();
        m_exposure = stream->readFloat();
        m_offset = stream->readFloat();
    }

    void configure() {
        if (m_nested == NULL)
            Log(EError, "Adjust plugin needs a nested texture!");
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture)))
            m_nested = static_cast<Texture *>(child);
        else
            Texture::addChild(name, child);
    }
    
    Spectrum eval(const Intersection &its, bool filter) const {
        Spectrum result = m_nested->eval(its, filter);
        std::string hashstr;
        if(its.instance != NULL) {
            hashstr = its.instance->getID();
        } else {
            hashstr = its.shape->getName();
        }
        Float value;
        std::map<std::string, float>* threadMap = g_threadMap.get();
        if( !threadMap ) {
            threadMap = new std::map<std::string, float>();
            g_threadMap.reset(threadMap);
        }
        auto it = threadMap->find(hashstr);
        if(it!=threadMap->end()) {
            value = threadMap->find(hashstr)->second;
        } else {
            ref<Random> random = new Random(m_hasher(hashstr));
            value = random->nextFloat() - 0.5f;
            threadMap->insert(std::pair<std::string, float>(hashstr,value));
        }
        // ref<Random> random = new Random(m_hasher(hashstr));
        // value = random->nextFloat() - 0.5f;

        if(m_gamma != 0.0f) {
            for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
                result[i] = pow(result[i], 1.0f/(1 + m_gamma * value));
            }
        }
        if(m_hue != 0.0f || m_saturation != 0.0f) {
            HSV hsv(result);
            hsv.h += m_hue * value;
            hsv.s += m_saturation * value;
            result = hsv.toSpectrum();
        }
        if(m_exposure != 0.0f) {
            for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
                result[i] *= pow(2.0f, m_exposure * value);
            }
        }
        if(m_offset != 0.0f) {
            result += Spectrum(m_offset * value);
        }
        return result;
    }

    void evalGradient(const Intersection &its, Spectrum *gradient) const {
        m_nested->evalGradient(its, gradient);
        gradient[0] *= m_scale;
        gradient[1] *= m_scale;
    }

    Spectrum getAverage() const {
        return m_nested->getAverage() * m_scale;
    }

    Spectrum getMaximum() const {
        return m_nested->getMaximum() * m_scale;
    }

    Spectrum getMinimum() const {
        return m_nested->getMinimum() * m_scale;
    }

    bool isConstant() const {
        return m_nested->isConstant();
    }

    ref<Bitmap> getBitmap(const Vector2i &sizeHint) const {
        ref<Bitmap> result = m_nested->getBitmap(sizeHint);

        if (m_scale == Spectrum(m_scale[0])) {
            result->scale(m_scale[0]);
        } else {
            result = result->convert(Bitmap::ESpectrum, Bitmap::EFloat);

            Spectrum *data = (Spectrum *) result->getFloatData();
            size_t pixelCount = result->getPixelCount();
            for (size_t i=0; i<pixelCount; ++i)
                *data++ *= m_scale;
        }

        return result;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "VariationTexture[" << endl
            << "  nested = " << indent(m_nested->toString()) << "," << endl
            << "  scale = " << m_scale.toString() << endl
            << "  gamma = " << m_gamma << endl
            << "  hue = " << m_hue << endl
            << "  saturation = " << m_saturation << endl
            << "  exposure = " << m_exposure << endl
            << "  offset = " << m_offset << endl
            << "]";
        return oss.str();
    }

    bool usesRayDifferentials() const {
        return m_nested->usesRayDifferentials();
    }

    bool isMonochromatic() const {
        return m_nested->isMonochromatic();
    }

    Shader *createShader(Renderer *renderer) const;

    void serialize(Stream *stream, InstanceManager *manager) const {
        Texture::serialize(stream, manager);
        manager->serialize(stream, m_nested.get());
        m_scale.serialize(stream);
        stream->writeFloat(m_gamma);
        stream->writeFloat(m_hue);
        stream->writeFloat(m_saturation);
        stream->writeFloat(m_exposure);
        stream->writeFloat(m_offset);
    }

    MTS_DECLARE_CLASS()
protected:
    ref<const Texture> m_nested;
    Spectrum m_scale;
    Float m_gamma;
    Float m_hue;
    Float m_saturation;
    Float m_exposure;
    Float m_offset;
    boost::hash<std::string> m_hasher;
};

// ================ Hardware shader implementation ================

class VariationTextureShader : public Shader {
public:
    VariationTextureShader(Renderer *renderer, const Texture *nested, const Spectrum &scale)
        : Shader(renderer, ETextureShader), m_nested(nested), m_scale(scale) {
        m_nestedShader = renderer->registerShaderForResource(m_nested.get());
    }

    bool isComplete() const {
        return m_nestedShader.get() != NULL;
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_nested.get());
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_nestedShader.get());
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_scale;" << endl
            << endl
            << "vec3 " << evalName << "(vec2 uv) {" << endl
            << "    return " << depNames[0] << "(uv) * " << evalName << "_scale;" << endl
            << "}" << endl;
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_scale", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &nestedUnitOffset) const {
        program->setParameter(parameterIDs[0], m_scale);
    }

    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_nested;
    ref<Shader> m_nestedShader;
    Spectrum m_scale;
};

Shader *VariationTexture::createShader(Renderer *renderer) const {
    return new VariationTextureShader(renderer, m_nested.get(), m_scale);
}

MTS_IMPLEMENT_CLASS(VariationTextureShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(VariationTexture, false, Texture2D)
MTS_EXPORT_PLUGIN(VariationTexture, "Adjust texture");
MTS_NAMESPACE_END
