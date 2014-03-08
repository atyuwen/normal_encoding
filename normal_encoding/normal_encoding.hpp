#pragma once

#include <cstdint>

struct float3
{
	float x, y, z;
};

// encode normalized 3d vector to 16-bits unsigned integer
uint16_t encode(const float3& v);

float3 decode(uint16_t s);
