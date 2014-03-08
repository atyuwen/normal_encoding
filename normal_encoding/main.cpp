#include <cstdlib>
#include <random>
#include <iostream>
#include <ctime>

#include "normal_encoding.hpp"

float clamp(float v, float low, float high)
{
	return v < low ? low : (v > high ? high : v);
}

int main()
{
#define PRECALC 1
#if PRECALC
	// calculate the parameters and look-up tables needed in encoding/decoding
	// the result is written to "out.txt".
	void precalc();
	precalc();
#endif


#define JAGGED_TO_LINEAR_TEST 0
#if JAGGED_TO_LINEAR_TEST
	void test();
	test();
#endif

	std::cout << std::endl;
	std::cout << "now testing..." << std::endl;
	srand(static_cast<unsigned int>(time(0)));

	const int num_test = 1000000;
	float max_err = 0;
	float avg_err = 0;
	for (int i = 0; i < num_test; ++i)
	{
		float x = static_cast<float>(rand()) / RAND_MAX - 0.5f;
		float y = static_cast<float>(rand()) / RAND_MAX - 0.5f;
		float z = static_cast<float>(rand()) / RAND_MAX - 0.5f;
		float len = sqrt(x * x + y * y + z * z);
		float3 v = {x / len, y / len, z / len};
		
		uint16_t s = encode(v);
		float3 u = decode(s);
		float dot = v.x * u.x + v.y * u.y + v.z * u.z;
		float err = acos(clamp(dot, -1.0f, 1.0f));
		if (err > max_err) max_err = err;
		avg_err += err;
	}
	avg_err /= num_test;

	const float c_pi = 3.14159265358979323846f;
	std::cout << "max error (in rad): " << max_err << std::endl;
	std::cout << "max error (in deg): " << max_err / c_pi * 180 << std::endl;
	std::cout << "avg error (in rad): " << avg_err << std::endl;
	std::cout << "avg error (in deg): " << avg_err / c_pi * 180 << std::endl;

	system("pause");
	return 0;
}
