#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <algorithm>

static std::default_random_engine engine(10);
static std::normal_distribution<double> normal(0, 1);

void sliced_matching(const double *source, double *target, int W, int H, double *image_result) {
	memcpy(image_result, source, W*H * 3 * sizeof(double));
	for(int t = 0; t < 100; ++t) {
		double xDir = normal(engine);
		double yDir = normal(engine);
		double zDir = normal(engine);
		double norm = sqrt(xDir * xDir + yDir * yDir + zDir * zDir);
		xDir /= norm;
		yDir /= norm;
		zDir /= norm;
		std::vector<std::pair<double, int>> result_sorted(W * H);
		std::vector<double> sortedTarget(W * H);
		for(int i = 0; i < W * H; ++i) {
			// Project image_result on target 
			double dotResult = image_result[i * 3] * xDir + image_result[i * 3 + 1] * yDir + image_result[i * 3 + 2] * zDir;
			double dotTarget = target[i * 3] * xDir + target[i * 3 + 1] * yDir + target[i * 3 + 2] * zDir;
			result_sorted[i] = std::make_pair(dotResult, i);
			sortedTarget[i] = dotTarget;
		}
		std::sort(result_sorted.begin(), result_sorted.end());
		std::sort(sortedTarget.begin(), sortedTarget.end());

		for(int i = 0; i < W * H; ++i) {
			double motionAmount = sortedTarget[i] - result_sorted[i].first;
			int index = result_sorted[i].second;
			image_result[index * 3] += motionAmount * xDir;
			image_result[index * 3 + 1] += motionAmount * yDir;
			image_result[index * 3 + 2] += motionAmount * zDir;
		}
	}
}

unsigned char clamp(double x, double min, double max) {
	if(x < min) return min;
	if(x > max) return max;
	return x;
}

int main() {

	int W, H, C;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *imageSouce = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);
	unsigned char *colorSource = stbi_load("redim.jpg",
								 &W,
								 &H,
								 &C,
								 STBI_rgb);
							 
	std::vector<double> imageSource(W*H*3);
	std::vector<double> imageTarget(W * H  * 3);
	for (int i=0; i<W*H*3; i++) {
		imageSource[i] = imageSouce[i];
		imageTarget[i] = colorSource[i];
	}

	std::vector<double> result(W*H*3);
	sliced_matching(&imageSource[0], &imageTarget[0], W, H, &result[0]);
	
	std::vector<unsigned char> image_result(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			image_result[(i*W + j) * 3] = clamp(result[(i*W + j) * 3], 0, 255);
			image_result[(i*W + j) * 3 + 1] = clamp(result[(i*W + j) * 3 + 1], 0, 255);
			image_result[(i*W + j) * 3 + 2] = clamp(result[(i*W + j) * 3 + 2], 0, 255);
		}
	}
	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

	return 0;
}