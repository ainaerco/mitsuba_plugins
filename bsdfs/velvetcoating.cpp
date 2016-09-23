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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/*!\plugin{ward}{Anisotropic Ward BRDF}
 * \order{15}
 * \parameters{
 *     \parameter{variant}{\String}{
 *         Determines the variant of the Ward model to use:
 *         \begin{enumerate}[(i)]
 *             \item \code{ward}: The original model by Ward \cite{Ward1992Measuring}
 *             --- suffers from energy loss at grazing angles.
 *             \item \code{ward-duer}: Corrected Ward model with lower energy loss
 *             at grazing angles \cite{Dur2006Improved}.
 *             Does not always conserve energy.
 *             \item \code{balanced}: Improved version of the \code{ward-duer}
 *             model with energy balance at all angles \cite{Geisler2010New}.
 *         \end{enumerate}
 *         Default: \texttt{balanced}
 *     }
 *     \parameter{alphaU, alphaV}{\Float\Or\Texture}{
 *         Specifies the anisotropic roughness values along the tangent and
 *         bitangent directions.
 *         \default{0.1}.
 *     }
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{
 *         Specifies the weight of the specular reflectance component.\default{0.2}}
 *     \parameter{diffuse\showbreak Reflectance}{\Spectrum\Or\Texture}{
 *         Specifies the weight of the diffuse reflectance component\default{0.5}}
 * }
 * \renderings{
 *     \rendering{$\alpha_u=0.1,\ \alpha_v=0.3$}{bsdf_ward_01_03}
 *     \rendering{$\alpha_u=0.3,\ \alpha_v=0.1$}{bsdf_ward_03_01}
 * }

 * This plugin implements the anisotropic Ward reflectance model and
 * several extensions. They are described in the papers
 * \begin{enumerate}[(i)]
 *    \item ``Measuring and Modeling Anisotropic Reflection''
 *      by Greg Ward \cite{Ward1992Measuring}
 *    \item ``Notes on the Ward BRDF'' by Bruce Walter \cite{Walter2005Notes}
 *    \item ``An Improved Normalization for the Ward Reflectance Model''
 *      by Arne D\"ur \cite{Dur2006Improved}
 *    \item ``A New Ward BRDF Model with Bounded Albedo'' by
 *      Geisler-Moroder et al. \cite{Geisler2010New}
 * \end{enumerate}
 *
 * Like the Phong BRDF, the Ward model does not take the Fresnel reflectance
 * of the material into account. In an experimental study by Ngan et al.
 * \cite{Ngan2005Experimental}, the Ward model performed noticeably worse than
 * models based on microfacet theory.
 *
 * For this reason, it is usually preferable to switch to a microfacet model
 * that incorporates knowledge about the material's index of refraction. In Mitsuba,
 * two such alternatives to \pluginref{ward} are given by the plugins
 * \pluginref{roughconductor} and \pluginref{roughplastic} (depending on the
 * material type).
 *
 * When using this plugin, note that the diffuse and specular reflectance
 * components should add up to a value less than or equal to one (for each
 * color channel). Otherwise, they will automatically be scaled appropriately
 * to ensure energy conservation.
 */
class VelvetCoating : public BSDF {
public:
	VelvetCoating(const Properties &props)
		: BSDF(props) {
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(0.2f)));

		Float alpha = props.getFloat("alpha", 0.1f),
			  alphaU = props.getFloat("alphaU", alpha),
			  alphaV = props.getFloat("alphaV", alpha);

		m_alphaU = new ConstantFloatTexture(alphaU);
		if (alphaU == alphaV)
			m_alphaV = m_alphaU;
		else
			m_alphaV = new ConstantFloatTexture(alphaV);
		m_specularSamplingWeight = 1.0f;
	}

	VelvetCoating(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_nested = static_cast<BSDF *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		configure();
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			if (m_nested != NULL)
				Log(EError, "Only a single nested BRDF can be added!");
			m_nested = static_cast<BSDF *>(child);
		} 
		else if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "alphaU")
				m_alphaU = static_cast<Texture *>(child);
			else if (name == "alphaV")
				m_alphaV = static_cast<Texture *>(child);
			else if (name == "specularReflectance")
				m_specularReflectance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
		manager->serialize(stream, m_nested.get());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
	}

	void configure() {
		if (!m_nested)
			Log(EError, "A child BSDF instance is required");
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV)
			extraFlags |= EAnisotropic;

		m_components.clear();
		for (int i=0; i<m_nested->getComponentCount(); ++i)
			m_components.push_back(m_nested->getType(i) | extraFlags);
		
		m_components.push_back(EGlossyReflection | EFrontSide | extraFlags
			| ((!m_specularReflectance->isConstant() || !m_alphaU->isConstant()
			  || !m_alphaV->isConstant()) ? ESpatiallyVarying : 0));

		m_usesRayDifferentials = m_nested->usesRayDifferentials()
			|| m_specularReflectance->usesRayDifferentials();
		BSDF::configure();
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		Spectrum result(0.0f);
		if (Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
			return result;
		bool hasSpecular = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool hasNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);

		Float probNested, probSpecular;
		if (hasSpecular && hasNested) {
			/* Find the probability of sampling the specular component */
			probSpecular = getProbability();
			probNested = 1 - probSpecular;
		} else {
			probNested = probSpecular = 1.0f;
		}

		if (hasSpecular) 
		{
			Vector H = normalize(bRec.wi+bRec.wo);
			Float alphaU = m_alphaU->eval(bRec.its).average();
			Float alphaV = m_alphaV->eval(bRec.its).average();
			const Float LdN = Frame::cosTheta(bRec.wo);
			const Float NdH = Frame::cosTheta(H);
			const Float VdH = dot(H, bRec.wi);
			const Float LdH = dot(H, bRec.wo);
			if(!(LdN<Epsilon || VdH<Epsilon || 
				NdH>1-Epsilon)) {
				const Float theta = std::acos(NdH);
				const Float sintheta = std::sinf(theta);

				const Float cottheta = NdH/sintheta;
				Float specRef = 1 / (4 * M_PI * alphaU * alphaV * 
					std::sqrt (VdH * LdH)) * std::exp(-cottheta*cottheta /(alphaU*alphaV));
				/* Important to prevent numeric issues when evaluating the
				   sampling density of the Ward model in places where it takes
				   on miniscule values (Veach-MLT does this for instance) */
				if (specRef > Epsilon)
					result += probSpecular * m_specularReflectance->eval(bRec.its) * specRef;
			}
		}
		if (hasNested) 
		{
			result += probNested * m_nested->eval(bRec, measure);
		}
		return result;
	}

	Float getProbability() const {
		return 0.5;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
			return 0.0f;

		bool sampleSpecular = (bRec.typeMask & EGlossyReflection)
			&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);

		Float probNested, probSpecular;
		if (sampleSpecular && sampleNested) {
			/* Find the probability of sampling the specular component */
			probSpecular = getProbability();
			probNested = 1 - probSpecular;
		} else {
			probNested = probSpecular = 1.0f;
		}

		Float pdf = 0.0f;
		if (sampleSpecular && Frame::cosTheta(bRec.wo) * Frame::cosTheta(bRec.wi) > 0) {
			Float alphaU = m_alphaU->eval(bRec.its).average();
			Float alphaV = m_alphaV->eval(bRec.its).average();
			Vector H = normalize(bRec.wi+bRec.wo);
			Float LdN = Frame::cosTheta(bRec.wo);
			Float NdH = Frame::cosTheta(H);
			Float VdH = dot(H, bRec.wi);
			if(!(LdN<Epsilon || VdH<Epsilon ||
				NdH>1-Epsilon)) {
				// Float cottheta = NdH / Frame::sinTheta(H);
				const Float theta = std::acos(NdH);
				const Float sintheta = std::sinf(theta);

				const Float cottheta = NdH/sintheta;
				pdf = cottheta/(4.0f * M_PI * alphaU * alphaU * VdH *
					std::powf(sintheta, 3)) * std::exp(-cottheta*cottheta /(alphaU*alphaU));
				pdf *= probSpecular;
				if (pdf < Epsilon)
					pdf = 0;
			}

		}
		if (sampleNested) {
			pdf += m_nested->pdf(bRec, measure) * probNested;
		}

		return pdf;
	}

	inline Spectrum sample(BSDFSamplingRecord &bRec, Float &_pdf, const Point2 &_sample) const {

		bool sampleSpecular = (bRec.typeMask & EGlossyReflection)
			&& (bRec.component == -1 || bRec.component == (int) m_components.size()-1);
		bool sampleNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
			&& (bRec.component == -1 || bRec.component < (int) m_components.size()-1);

		if (!sampleSpecular && !sampleNested)
			return Spectrum(0.0f);
		
		bool choseSpecular = sampleSpecular;
		Float probSpecular = getProbability();

		Point2 sample(_sample);
		if (sampleSpecular && sampleNested) {
			if (sample.x < probSpecular) {
				sample.x /= probSpecular;
			} else {
				sample.x = (sample.x - probSpecular) / (1 - probSpecular);
				choseSpecular = false;
			}
		}
		if (choseSpecular) {
			Float alphaU = m_alphaU->eval(bRec.its).average();
			Float alphaV = m_alphaV->eval(bRec.its).average();
			Float phiH = 2.0f * M_PI * sample.y;
			Float thetaH = std::atan(1/(alphaU)*math::safe_sqrt(-math::fastlog(sample.x)));
			Vector H = sphericalDirection(thetaH, phiH);
			bRec.wo = H * (2.0f * dot(bRec.wi, H)) - bRec.wi;
			bRec.sampledComponent = (int) m_components.size() - 1;
			bRec.sampledType = EGlossyReflection;
			bRec.eta = 1.0f;
			if (Frame::cosTheta(bRec.wo) <= 0.0f)
				return Spectrum(0.0f);
		} else {
			Spectrum result = m_nested->sample(bRec, _pdf, sample);
			if (result.isZero())
				return Spectrum(0.0f);

		}
		EMeasure measure = getMeasure(bRec.sampledType);
		_pdf = pdf(bRec, measure);
		if (_pdf == 0)
			return Spectrum(0.0f);
		else
			return eval(bRec, measure) / _pdf;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		Float pdf;
		return VelvetCoating::sample(bRec, pdf, sample);
	}

	Float getRoughness(const Intersection &its, int component) const {
		return component < (int) m_components.size() - 1
			? m_nested->getRoughness(its, component)
			: 0.5f * (m_alphaU->eval(its).average()
				+ m_alphaV->eval(its).average());
	}

	Shader *createShader(Renderer *renderer) const;

	std::string toString() const {
		std::ostringstream oss;
		oss << "VelvetCoating[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
			<< "  alphaV = " << indent(m_alphaV->toString()) << endl
			<< "  nested = " << indent(m_nested.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU;
	ref<Texture> m_alphaV;
	ref<BSDF> m_nested;
	Float m_specularSamplingWeight;
};

// ================ Hardware shader implementation ================

/**
 * GLSL port of the Ward shader. This version only implements the variant
 * with energy balance. When the roughness is lower than
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview.
 */
class VelvetCoatingShader : public Shader {
public:
	VelvetCoatingShader(Renderer *renderer, const BSDF *nested,
			const Texture *specularColor,
			const Texture *alphaU,
			const Texture *alphaV) : Shader(renderer, EBSDFShader),
			m_nested(nested),
			m_specularReflectance(specularColor),
			m_alphaU(alphaU), m_alphaV(alphaV) {
		m_nestedShader = renderer->registerShaderForResource(m_nested.get());
		m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
		m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
		m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());
	}

	bool isComplete() const {
		return m_nestedShader.get() != NULL &&
			   m_specularReflectanceShader.get() != NULL &&
			   m_alphaU.get() != NULL &&
			   m_alphaV.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_nestedShader.get());
		deps.push_back(m_specularReflectanceShader.get());
		deps.push_back(m_alphaUShader.get());
		deps.push_back(m_alphaVShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_nested.get());
		renderer->unregisterShaderForResource(m_specularReflectance.get());
		renderer->unregisterShaderForResource(m_alphaU.get());
		renderer->unregisterShaderForResource(m_alphaV.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z <= 0.0 || wo.z <= 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    vec3 H = normalize(wi + wo);" << endl
			<< "    float alphaU = max(0.3, " << depNames[2] << "(uv)[0]);" << endl
			<< "    float alphaV = max(0.3, " << depNames[3] << "(uv)[0]);" << endl
			<< "	float cottheta = cosTheta(H)/sinTheta(H);" << endl
			<< "    float factor1 = 1/(4*3.1415*alphaU*alphaV*sqrt(dot(H, wi) * dot(H, wo) ));"  << endl
			<< "    float exponent = -cottheta*cottheta /(alphaU*alphaV);" << endl
			<< "    float specRef = factor1 * exp(exponent);" << endl
			<< "    vec3 nested = " << depNames[0] << "(uv, wi, wo);" << endl
			<< "    return " << depNames[1] << "(uv) * specRef + nested;" << endl
			<< "}" << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << depNames[0] << "_diffuse(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const BSDF> m_nested;
	ref<const Texture> m_specularReflectance;
	ref<const Texture> m_alphaU;
	ref<const Texture> m_alphaV;
	
	ref<Shader> m_nestedShader;
	ref<Shader> m_specularReflectanceShader;
	ref<Shader> m_alphaUShader;
	ref<Shader> m_alphaVShader;
};

Shader *VelvetCoating::createShader(Renderer *renderer) const {
	return new VelvetCoatingShader(renderer, m_nested.get(),
		m_specularReflectance.get(), m_alphaU.get(), m_alphaV.get());
}

MTS_IMPLEMENT_CLASS(VelvetCoatingShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(VelvetCoating, false, BSDF);
MTS_EXPORT_PLUGIN(VelvetCoating, "Anisotropic Ward BRDF");
MTS_NAMESPACE_END
