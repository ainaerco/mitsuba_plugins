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
#include <sstream>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/core/fresolver.h>
#include <OpenImageIO/texture.h>
#include "OpenImageIO/sysutil.h"
#include "OpenImageIO/strutil.h"

#include <boost/algorithm/string.hpp>

using namespace OpenImageIO;

MTS_NAMESPACE_BEGIN

class OiioTexture : public Texture2D {
public:
	typedef TSpectrum<Float, 3> Color3;
	~OiioTexture()
	{
		Log(EInfo, "Texture memory usage: %s", 
			Strutil::memformat(Sysutil::memory_used(true)));
		TextureSystem::destroy (m_texsys);
	}

	OiioTexture(const Properties &props) : Texture2D(props) {
		fs::path filePath = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		// if (!fs::exists(filePath))
		// 	Log(EError, "Texture file \"%s\" could not be found!", filePath.string().c_str());
		m_filename = filePath.string();
		// m_ustring = ustring(filePath.string().c_str());
		m_texsys = TextureSystem::create();
		// m_texsys->attribute ("statistics:level", 2);
		// m_texsys->attribute ("autotile", 1);
		m_texsys->attribute ("automip", 1);

		// m_perthread_info = m_texsys->get_perthread_info();
		// m_texture_handle = m_texsys->get_texture_handle(m_ustring);
		m_defaultColor = props.getSpectrum("defaultColor", Spectrum(.5f));
		// m_spec = m_texsys->imagespec(m_ustring);
		float sblur = 0;
		float width = 1;
		m_options.sblur = sblur;
		m_options.tblur = sblur;
		m_options.swidth = width;
		m_options.twidth = width;
		m_options.interpmode = TextureOpt::InterpSmartBicubic;
		m_options.mipmode = TextureOpt::MipModeDefault;
		m_options.twrap = m_options.swrap = TextureOpt::WrapPeriodic;
		m_options.firstchannel = 0;
		
	}

	OiioTexture(Stream *stream, InstanceManager *manager)
	 : Texture2D(stream, manager) {
		m_filename = stream->readString();
		Log(EDebug, "Unserializing oiio texture \"%s\"", m_filename.c_str());
		m_ustring = ustring(m_filename.c_str());
		m_texsys = TextureSystem::create();
		m_spec = m_texsys->imagespec(m_ustring);
		m_defaultColor = Spectrum(stream);
		m_options.sblur = stream->readFloat();
		m_options.tblur = stream->readFloat();
		m_options.swidth = stream->readFloat();
		m_options.twidth = stream->readFloat();
		m_options.interpmode = (TextureOpt::InterpMode)stream->readUInt();
		m_options.mipmode = (TextureOpt::MipMode)stream->readUInt();
		m_options.swrap = (TextureOpt::Wrap)stream->readUInt();
		m_options.twrap = (TextureOpt::Wrap)stream->readUInt();
		m_options.firstchannel = stream->readUInt();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture2D::serialize(stream, manager);
		stream->writeString(m_filename);
		m_defaultColor.serialize(stream);
		stream->writeFloat(m_options.sblur);
		stream->writeFloat(m_options.tblur);
		stream->writeFloat(m_options.swidth);
		stream->writeFloat(m_options.twidth);
		stream->writeUInt(m_options.interpmode);
		stream->writeUInt(m_options.mipmode);
		stream->writeUInt(m_options.swrap);
		stream->writeUInt(m_options.twrap);
		stream->writeUInt(m_options.firstchannel);
	}

	inline Spectrum eval(const Point2 &uv) const {
		std::stringstream ss;
		ss << 1001 + 10*int(-uv.y) + int(uv.x);
		std::string filename = m_filename;
		boost::replace_last(filename, "<UDIM>", ss.str());
		ustring ustr = ustring(filename.c_str());
		const ImageSpec* spec = m_texsys->imagespec(ustr);
		float *sampleResult = ALLOCA (float, spec->nchannels);
		float dsdx = (1.0f / spec->width) * 1.5f;
		float dtdy = (1.0f / spec->height) * 1.5f;
		bool ok = m_texsys->texture (ustr, const_cast<TextureOpt&>(m_options), 
			uv.x, uv.y, dsdx, 0.0f, 0.0f, dtdy, 
			spec->nchannels, sampleResult, NULL, NULL);
		Spectrum result;
		if(!ok)
			return result;
		result.fromLinearRGB(sampleResult[0], sampleResult[1], sampleResult[2]);
		return result;

	}

	Spectrum eval(const Point2 &uv,
			const Vector2 &d0, const Vector2 &d1) const {
		std::stringstream ss;
		ss << 1001 + 10*int(-uv.y) + int(uv.x);
		std::string filename = m_filename;
		boost::replace_last(filename, "<UDIM>", ss.str());
		ustring ustr = ustring(filename.c_str());
		const ImageSpec* spec = m_texsys->imagespec(ustr);
		float *sampleResult = ALLOCA (float, spec->nchannels);
		bool ok = m_texsys->texture (ustr, const_cast<TextureOpt&>(m_options),
			uv.x, uv.y, d0.x/100.0f, d0.y/100.0f, d1.x/100.0f, d1.y/100.0f, 
			spec->nchannels, sampleResult, NULL, NULL);
		Spectrum result;
		if(!ok)
			return result;
		result.fromLinearRGB(sampleResult[0], sampleResult[1], sampleResult[2]);
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
		float avg[4];
    	bool ok = m_texsys->get_texture_info (m_ustring, 0, ustring("averagecolor"),
                                   TypeDesc(TypeDesc::FLOAT,4), avg);
    	Spectrum result;
    	if(ok)
    		result.fromLinearRGB(avg[0],avg[1],avg[2]);
		return result;
	}

	bool isConstant() const {
		return false;
		float avg[4];
		return m_texsys->get_texture_info (m_ustring, 0, ustring("constantcolor"),
                                   TypeDesc(TypeDesc::FLOAT,4), avg);
	}

	bool isMonochromatic() const {
		return false;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "OiioTexture[" << endl
			<< "    filename = " << m_filename << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	std::string m_filename;
	ustring m_ustring;
	TextureSystem *m_texsys;
	TextureSystem::TextureHandle *m_texture_handle;
	TextureSystem::Perthread *m_perthread_info;
	TextureOpt m_options;
	Spectrum m_defaultColor;
	const ImageSpec* m_spec;
};

// ================ Hardware shader implementation ================

class OiioTextureShader : public Shader {
public:
	OiioTextureShader(Renderer *renderer, const Spectrum &value)
		: Shader(renderer, ETextureShader), m_value(value) {
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_value;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return " << evalName << "_value;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_value", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_value);
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_value;
};

Shader *OiioTexture::createShader(Renderer *renderer) const {
	return new OiioTextureShader(renderer, m_defaultColor);
}

MTS_IMPLEMENT_CLASS(OiioTextureShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(OiioTexture, false, Texture2D)
MTS_EXPORT_PLUGIN(OiioTexture, "OiioTexture texture");
MTS_NAMESPACE_END
