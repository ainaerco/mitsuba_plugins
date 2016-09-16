#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/platform.h>

MTS_NAMESPACE_BEGIN

class HSV
{
public:
    HSV(){}
    HSV(Spectrum c)
    {
        Float r,g,b;
        c.toLinearRGB(r,g,b);
        Float cMin = std::min(r, std::min(g, b));
        Float cMax = std::max(r, std::max(g, b));
        Float delta = cMax - cMin; 
        v = cMax;
        if (cMax > 0)
            s = delta / cMax;
        else s = 0;
        if (s <= 0)
            h = 0;
        else
        {
            if      (r >= cMax) h = (g-b) / delta;
            else if (g >= cMax) h = 2 + (b-r) / delta;
            else                h = 4 + (r-g) / delta;
            h /= 6;
            h += h < 0;
        }
    }

    inline Spectrum toSpectrum()
    {
        Spectrum result;
        Float r,g,b;
        if (s < 0.0001f)
        {
            r = g = b = v;
        }
        else
        {
            h = 6.0f * (h - floorf(h));
            int hi = (int) h;
            float f = h - hi;
            float p = v * (1.0f - s);
            float q = v * (1.0f - s * f);
            float t = v * (1.0f - s * (1.0f - f));
            switch (hi)
            {
                case 0 : r = v; g = t; b = p;
                break;
                case 1 : r = q; g = v; b = p;
                break;
                case 2 : r = p; g = v; b = t;
                break;
                case 3 : r = p; g = q; b = v;
                break;
                case 4 : r = t; g = p; b = v;
                break;
                default: r = v; g = p; b = q;
                break;
            }
        }
        result.fromLinearRGB(r,g,b);
        return result;
    }

    Float h,s,v;
};

MTS_NAMESPACE_END