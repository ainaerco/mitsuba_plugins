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
#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/basicshader.h>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

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

class BitmapUdim : public Texture {
public:
	BitmapUdim(const Properties &props) : Texture(props), m_count(0) {
		std::cout<<"BitmapUdim props\n";


		if (!props.hasProperty("filename"))
			Log(EError, "'filename' parameter shoud be specified!");
		m_filename = props.getString("filename");
		std::size_t found = m_filename.rfind("####");

		if (found == std::string::npos)
			Log(EError, "Use #### in filename to use udims!");

		boost::filesystem::path filepath(m_filename);
		boost::filesystem::path extension = filepath.extension();
		boost::filesystem::path dir = filepath.parent_path();
		boost::filesystem::directory_iterator b(dir), e;
		boost::regex reg("\\.\\d\\d\\d\\d");

		for (auto i=b; i!=e; ++i)
		{
			const boost::filesystem::path& path = i->path();
			std::string udimStr = path.stem().extension().string();
			if(extension!=path.extension() || 
				m_filename.length()!=path.string().length() ||
				path.string().substr(0,found-1)!=m_filename.substr(0,found-1) ||
				udimStr.length()!=5)
				continue;
			// if(boost::regex_search(path.stem().string().substr(found,path.stem().string().length()), reg))
			if (!boost::regex_search(udimStr, reg))
				continue;
			int udim = boost::lexical_cast<int>( udimStr.substr(1,udimStr.length()-1) );
			Properties bitmapProps("bitmap");
			bitmapProps.setPluginName("bitmap");
			bitmapProps.setString("filename", path.string());
			Texture* bitmap = static_cast<Texture *> (PluginManager::getInstance()->createObject(MTS_CLASS(Texture), bitmapProps));
			m_udims.push_back(udim);
			m_bitmaps.push_back(bitmap);
			m_count++;
		}
		if (m_count==0)
			Log(EError, "Textures count is 0!");
	}

	BitmapUdim(Stream *stream, InstanceManager *manager)
		: Texture(stream, manager) {
		m_count = stream->readSize();
		for (size_t i=0; i<m_count; ++i) {
			m_udims.push_back(stream->readInt());
			Texture *texture = static_cast<Texture *>(manager->getInstance(stream));
			texture->incRef();
			m_bitmaps.push_back(texture);
		}

	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			Texture *texture = static_cast<Texture *>(child);
			m_bitmaps.push_back(texture);
			texture->incRef();
		} else {
			Texture::addChild(name, child);
		}
	}
	
	Spectrum eval(const Intersection &its, bool filter) const {
		int udim = 1001 + 10*int(-its.uv.y) + int(its.uv.x);
		auto it = std::find(m_udims.begin(),m_udims.end(),udim);
		if(it==m_udims.end()){
			return Spectrum();
		}
		return m_bitmaps[it - m_udims.begin()]->eval(its, filter);
	}

	void evalGradient(const Intersection &its, Spectrum *gradient) const {
		int udim = 1001 + 10*int(-its.uv.y) + int(its.uv.x);
		auto it = std::find(m_udims.begin(),m_udims.end(),udim);
		if(it==m_udims.end()){
			return ;
		}
		m_bitmaps[it - m_udims.begin()]->evalGradient(its, gradient);
	}

	Spectrum getAverage() const {
		Spectrum average(0.0f);
		for (size_t i=0; i<m_count; ++i) {
			average = average + m_bitmaps[i]->getAverage();
		}
		average = average/float(m_count);
		return average;
	}

	Spectrum getMaximum() const {
		Spectrum max = m_bitmaps[0]->getMaximum();
		for (size_t i=1; i<m_count; ++i) {
			Spectrum maxt = m_bitmaps[i]->getMaximum();
			for (int j=0; j<SPECTRUM_SAMPLES; ++j)
				max[j] = std::max(max[j], maxt[j]);
		}
		return max;
	}

	Spectrum getMinimum() const {
		Spectrum min = m_bitmaps[0]->getMinimum();
		for (size_t i=1; i<m_count; ++i) {
			Spectrum mint = m_bitmaps[i]->getMinimum();
			for (int j=0; j<SPECTRUM_SAMPLES; ++j)
				min[j] = std::min(min[j], mint[j]);
		}
		return min;
	}

	bool isConstant() const {
		bool isConstant = true;
		for (size_t i=0; i<m_count; ++i) {
			isConstant = isConstant && m_bitmaps[i]->isConstant();
		}
		return isConstant;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "BitmapUdim[" << endl
			<< "  filename = " << m_filename << "," << endl
			<< "]";
		return oss.str();
	}

	bool usesRayDifferentials() const {
		bool usesRayDifferentials = false;
		for (size_t i=0; i<m_count; ++i) {
			usesRayDifferentials = usesRayDifferentials || m_bitmaps[i]->usesRayDifferentials();
		}
		return usesRayDifferentials;
	}

	bool isMonochromatic() const {
		bool isMonochromatic = true;
		for (size_t i=0; i<m_count; ++i) {
			isMonochromatic = isMonochromatic && m_bitmaps[i]->isMonochromatic();
		}
		return isMonochromatic;
	}

	Shader *createShader(Renderer *renderer) const;

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		stream->writeSize(m_count);
		for (size_t i=0; i<m_count; ++i) {
			stream->writeInt(m_udims[i]);
			manager->serialize(stream, m_bitmaps[i].get());
		}
	}

	MTS_DECLARE_CLASS()
protected:
	std::string m_filename;
	std::vector<ref<Texture>> m_bitmaps;
	std::vector<int> m_udims;
	size_t m_count;
};

// ================ Hardware shader implementation ================

class BitmapUdimShader : public Shader {
public:
	BitmapUdimShader(Renderer *renderer, const Texture *nested)
		: Shader(renderer, ETextureShader), m_nested(nested) {
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
		oss << "uniform vec3 " << evalName << "_shader;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return " << depNames[0] << "(uv);" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &nestedUnitOffset) const {
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_nested;
	ref<Shader> m_nestedShader;
};

Shader *BitmapUdim::createShader(Renderer *renderer) const {
	return new BitmapUdimShader(renderer, m_bitmaps[0].get());
}

MTS_IMPLEMENT_CLASS(BitmapUdimShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(BitmapUdim, false, Texture2D)
MTS_EXPORT_PLUGIN(BitmapUdim, "Adjust texture");
MTS_NAMESPACE_END
