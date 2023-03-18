// Copyright © 2008-2022 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#pragma once

#include <vector>
#include "vector2.h"
#include <math.h>
#include <cassert>
#include <array>
#include "Random.h"
#include <iostream>

// distance between 2d points with an error of 5%
static double pseudo_dist(const vector2d &p1, const vector2d &p2)
{
	double d1 = abs(p1.x - p2.x);
	double d2 = abs(p1.y - p2.y);
	return 0.668 * (d1 + d2) + 0.332 * abs(d1 - d2);
}

static double tj_old(double ti, vector2d Pi, vector2d Pj){
	double t = sqrt(sqrt( pow((Pj.x - Pi.x), 2) + pow((Pj.y - Pi.y), 2)) ) + ti;
	return t;
}

static double tj(double ti, const vector2d &Pi, const vector2d &Pj)
{
	return sqrt(pseudo_dist(Pi, Pj)) + ti;
}

// based on
// https://www.programmersought.com/article/52274393478/

// between the point at the given index and the next, insert the given number of spline points
// the previous point is also used, and the point after the next
// we assume that the array is looped
static void catmull_rom_spline_4points(std::vector<vector2d> &spline, int prev, int index, int count)
{
	unsigned size = spline.size();

	const vector2d &P0 = spline[prev];
	const vector2d &P1 = spline[index];
	const vector2d &P2 = spline[(index + 1) % size];
	const vector2d &P3 = spline[(index + 2) % size];

	double t0 = 0;
	double t1 = tj(t0, P0, P1);
	double t2 = tj(t1, P1, P2);
	double t3 = tj(t2, P2, P3);

	// Can be understood as the interval between points
	double linespace = (t2 - t1) / (count + 1);

	if (linespace == 0) return;

	std::vector<vector2d> subSpline;
	subSpline.reserve(count);
	double t = t1 + linespace;
	for(int i = 0; i < count; ++i ){
		double A1_x = (t1-t)/(t1-t0)*P0.x + (t-t0)/(t1-t0)*P1.x;
		double A1_y = (t1-t)/(t1-t0)*P0.y + (t-t0)/(t1-t0)*P1.y;
		double A2_x = (t2-t)/(t2-t1)*P1.x + (t-t1)/(t2-t1)*P2.x;
		double A2_y = (t2-t)/(t2-t1)*P1.y + (t-t1)/(t2-t1)*P2.y;
		double A3_x = (t3-t)/(t3-t2)*P2.x + (t-t2)/(t3-t2)*P3.x;
		double A3_y = (t3-t)/(t3-t2)*P2.y + (t-t2)/(t3-t2)*P3.y;
		double B1_x = (t2-t)/(t2-t0)*A1_x + (t-t0)/(t2-t0)*A2_x;
		double B1_y = (t2-t)/(t2-t0)*A1_y + (t-t0)/(t2-t0)*A2_y;
		double B2_x = (t3-t)/(t3-t1)*A2_x + (t-t1)/(t3-t1)*A3_x;
		double B2_y = (t3-t)/(t3-t1)*A2_y + (t-t1)/(t3-t1)*A3_y;
		double C_x = (t2-t)/(t2-t1)*B1_x + (t-t1)/(t2-t1)*B2_x;
		double C_y = (t2-t)/(t2-t1)*B1_y + (t-t1)/(t2-t1)*B2_y;
		C_x = floor(C_x);
		C_y = floor(C_y);
		subSpline.push_back({ C_x, C_y });
		t = t + linespace;
	}
	spline.insert(spline.begin() + index + 1, subSpline.begin(), subSpline.end());
}

static void noise_between_points(std::vector<vector2d> &spline, int index, int count, Random &rand) {
	const float maxOffsetK = 0.6;
	std::vector<vector2d> noise;
	vector2d &p1 = spline[index];
	vector2d &p2 = spline[(index + 1) % spline.size()];
	vector2d along{ (p2.x - p1.x) / (count + 1), (p2.y - p1.y) / (count + 1) };
	// ortho vector: -b, a
	vector2d across{ -maxOffsetK * along.y, maxOffsetK * along.x };
	for(int i = 0; i < count; ++i) {
		vector2d newpnt = p1 + along * (i + 1) + across * rand.Double(-1.0, 1.0);
		noise.push_back(newpnt);
	}
	spline.insert(spline.begin() + index + 1, noise.begin(), noise.end());
}

template<bool put = true>
static void put_point(uint8_t *bitset, uint32_t x, uint32_t y, uint32_t rows, uint32_t cols)
{
	uint32_t bitset_idx = x >> 3;
	uint8_t bitmask = 1 << (x & 7);
	uint8_t &cell = *(bitset + bitset_idx + y * cols);
	if constexpr (put) {
		cell |= bitmask;
	} else {
		cell &= ~bitmask;
	}
}

static bool check_point(uint8_t *bitset, int x, int y, uint32_t rows, uint32_t cols)
{
	if (x < 0 || y < 0 || x >= (int)rows || y >= (int)rows) return true;
	uint32_t bitset_idx = x >> 3;
	uint8_t bitmask = 1 << (x & 7);
	uint8_t &cell = *(bitset + bitset_idx + y * cols);
	return (cell & bitmask) != 0;
}


static void line_between(uint8_t *bitset, vector2d &p1, vector2d &p2, uint32_t rows, uint32_t cols) {
	int dist = pseudo_dist(p1, p2);
	vector2d along{ (p2.x - p1.x) / (dist + 1), (p2.y - p1.y) / (dist + 1) };
	for(int i = 0; i < dist; ++i) {
		// -1 .. 1
		vector2d newpnt = p1 + along * (i + 1);
		put_point(bitset, newpnt.x, newpnt.y, rows, cols);
	}
}

static void closed_noisy_spline(std::vector<vector2d> &nodes, int count, float minsize, Random &rand) {
	bool updated = false;
	unsigned i, prev;
	do {
		i = 0;
		prev = nodes.size() - 1;
		updated = false;
		// spline stadia
		do {
			if (pseudo_dist(nodes[i], nodes[(i + 1) % nodes.size()]) > minsize) {
				catmull_rom_spline_4points(nodes, prev, i, count);
				updated = true;
			}
			prev = i;
			i += count + 1;
		} while (i < nodes.size());
		// noise stadia
		i = 0;
		do {
			if (pseudo_dist(nodes[i], nodes[(i + 1) % nodes.size()]) > minsize) {
				noise_between_points(nodes, i, count, rand);
				updated = true;
			}
			i += count + 1;
		} while (i < nodes.size());
	} while(updated);
}

static void closed_spline(std::vector<vector2d> &nodes, int count, float minsize, Random &rand) {
	bool updated = false;
	unsigned i, prev;
	do {
		i = 0;
		prev = nodes.size() - 1;
		updated = false;
		// spline stadia
		do {
			if (pseudo_dist(nodes[i], nodes[(i + 1) % nodes.size()]) > minsize) {
				catmull_rom_spline_4points(nodes, prev, i, count);
				updated = true;
			}
			prev = i;
			i += count + 1;
		} while (i < nodes.size());
	} while(updated);
}

static void fill(uint8_t *bitset, uint32_t rows, uint32_t cols, int x, int y)
{
	using point = std::array<int, 2>;
	std::vector<point> lifo;
	lifo.push_back({ x, y });
	put_point(bitset, x, y, rows, cols);
	while (!lifo.empty()) {
		auto curr = lifo.back();
		lifo.pop_back();
		point deltas[4] = {{-1, 0}, { 1, 0 }, { 0, -1}, { 0, 1 }};
		for(const auto delta : deltas) {
			auto x = curr[0] + delta[0];
			auto y = curr[1] + delta[1];
			if (!check_point(bitset, x, y, rows, cols)) {
				put_point(bitset, x, y, rows, cols);
				lifo.push_back({ x, y });
			}
		}
	}
}

struct Boolpnt {
	vector2d pos;
	bool white;
};

struct IndexDist {
	int index;
	double dist;
};

std::pair<IndexDist, IndexDist> two_nearest(std::vector<Boolpnt> &centers, const vector2d &pnt)
{
	IndexDist best { -1, std::numeric_limits<double>::max() };
	IndexDist after;

	// looking for the nearest by brute force
	for (int i = 0; i < centers.size(); ++i) {
		double dist = pseudo_dist(centers[i].pos, pnt);
		if (dist < best.dist) {
			after = best;
			best = { i, dist };
		}
	}

	return { best, after };
}

std::pair<IndexDist, IndexDist> two_nearest_from_indices(std::vector<Boolpnt> &centers, int near_indices[9], const vector2d &pnt)
{
	IndexDist best { -1, std::numeric_limits<double>::max() };
	IndexDist after;

	// looking for the nearest by brute force
	for (int i = 0; i < 9; ++i) {
		if (near_indices[i] == -1) continue;
		double dist = pseudo_dist(centers[near_indices[i]].pos, pnt);
		if (dist < best.dist) {
			after = best;
			best = { near_indices[i], dist };
		} else if (dist < after.dist) {
			after = { near_indices[i], dist };
		}
	}

	return { best, after };
}

static void create_patches(uint8_t *bitset, uint32_t rows, uint32_t cols, float size, float gap, Random &rand)
{
	int grid_cols = cols * 8 / size;
	int grid_rows = rows / size;
	// an odd number of grid elements is desirable so that the station is in the center of the square
	if (!(grid_cols & 1)) ++grid_cols;
	if (!(grid_rows & 1)) ++grid_rows;

	double x_step = cols * 8 / (double)grid_cols;
	double y_step = rows / (double)grid_rows;
	double margin = size * 0.1;
	// for a quick lookup, we will divide area into cells, in each cell there will be only one point
	std::vector<Boolpnt> centers;
	centers.resize(grid_cols * grid_rows);
	for (int y = 0; y < grid_rows; ++y) {
		for (int x = 0; x < grid_cols; ++x) {
			int i = y * grid_cols + x;
			double lx = rand.Double(margin, x_step - margin);
			double ly = rand.Double(margin, y_step - margin);
			centers[i].pos.x = x * x_step + lx;
			centers[i].pos.y = y * y_step + ly;
			centers[i].white = check_point(bitset, centers[i].pos.x, centers[i].pos.y, rows, cols);
		}
	}

	// random black patches
	int pcount = rand.Int32(4, 9);
	int maxtry = pcount * 5;

	for (int attempt = 0; attempt < maxtry; ++attempt) {
		int toblack = rand.Int32(0, centers.size() - 1);
		if (centers[toblack].white) {
			centers[toblack].white = false;
			--pcount;
			if (pcount <= 0) break;
		}
	}

	// the station must be on a white patch
	int station_i = grid_rows * (grid_cols / 2) + grid_rows / 2;
	centers[station_i].pos.x = grid_cols * x_step * 0.5;
	centers[station_i].pos.y = grid_rows * y_step * 0.5;
	centers[station_i].white = true;


	// clear bitset
	std::memset(bitset, 0, rows * cols);

	// stroke patches as voronoi diagram
	for (int y = 0; y < rows; ++y) {
		for (int x = 0; x < cols * 8; ++x) {
			int cell_x = std::floor(x / x_step);
			int cell_y = std::floor(y / y_step);
			int near_steps[9][2] = { {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {0, 0} };
			int near_indices[9];
			for (int i = 0; i < 9; ++i) {
				int near_x = cell_x + near_steps[i][0];
				int near_y = cell_y + near_steps[i][1];
				if (near_x >= 0 && near_x < grid_cols && near_y >= 0 && near_y < grid_rows) {
					near_indices[i] = near_y * grid_cols + near_x;
				} else {
					near_indices[i] = -1;
				}
			}

			auto nears = two_nearest_from_indices(centers, near_indices, vector2d(x, y));
			//auto nears = two_nearest(centers, vector2d(x, y));
			 if (centers[nears.first.index].white && std::abs(nears.first.dist - nears.second.dist) > gap) {
	//		if (std::abs(nears.first.dist - nears.second.dist) > gap) {
				put_point(bitset, x, y, rows, cols);
			}
		}
	}

		/*
	for (auto &p: centers) {
		put_point<false>(bitset, p.pos.x, p.pos.y, rows, cols);
	}

	for (int x = 0; x < grid_cols; ++x) {
		for (int y = 0; y < rows; ++y) {
			put_point<false>(bitset, x * x_step, y, rows, cols);
		}
	}

	for (int x = 0; x < cols * 8; ++x) {
		for (int y = 0; y < grid_rows; ++y) {
			put_point<false>(bitset, x, y * y_step, rows, cols);
		}
	}
	*/


}

static void create_patches_old(uint8_t *bitset, uint32_t rows, uint32_t cols, float size, float gap, Random &rand)
{
	// estimated amount
	float some_empirical_ratio = 0.3;
	int pcount = rows * cols * 8 / size / size * some_empirical_ratio;
	int maxtry = pcount * 8;

	// we will remember whether the point has already been painted over
	std::vector<Boolpnt> centers;

	int attempt;

	// spaceport position an center
	centers.push_back({ vector2d(cols * 8 / 2, rows / 2), true });
	// generate patch centers
	for (attempt = 0; attempt < maxtry; ++attempt) {
		vector2d newp(rand.Int32(gap * 0.5, cols * 8 - gap * 0.5), rand.Int32(gap * 0.5, rows - gap * 0.5));
		bool can_insert = true;
		for (auto &p : centers) {
			if (pseudo_dist(p.pos, newp) < size) {
				can_insert = false;
				break;
			}
		}
		if (can_insert) {
			bool white = check_point(bitset, newp.x, newp.y, rows, cols);
			centers.push_back({ newp, white });
			if (centers.size() >= pcount) break;
		}
	}

	printf("pcount: %d maxtry: %d attempts: %d inserted: %d effect: %.4f\n", pcount, maxtry, attempt, (int)centers.size(), (float)centers.size()/attempt);
	fflush(stdout);


	assert(centers.size() >= 2 && "we will not be able to search for the closest 2 if the array is less than 2");

	// random black patches
	pcount = rand.Int32(0, 3);
	maxtry = pcount * 5;

	for (attempt = 0; attempt < maxtry; ++attempt) {
		int toblack = rand.Int32(1, centers.size());
		if (centers[toblack].white) {
			centers[toblack].white = false;
			--pcount;
			if (pcount <= 0) break;
		}
	}

	// clear bitset
	std::memset(bitset, 0, rows * cols);

	// fill patches as voronoi diagram
	for (int x = 0; x < cols * 8; ++x) {
		for (int y = 0; y < rows; ++ y) {
			auto nears = two_nearest(centers, vector2d(x, y));
			if (centers[nears.first.index].white && std::abs(nears.first.dist - nears.second.dist) > gap) {
				put_point(bitset, x, y, rows, cols);
			}
		}
	}

}

// initial idea taken from
// https://plottersvg.ru/generator-spot
static void generate_blob(uint8_t *bitset, uint32_t seed, uint32_t rows, uint32_t cols, uint32_t points, uint32_t parts, double mindist)
{
	Random rand(seed);
	std::memset(bitset, 0, rows * cols);
	std::vector<vector2d> nodes;
	float max_radius = rows / 2;

	int cx = rows / 2;
	int cy = rows / 2;

	// random points around the center
	for (uint32_t i = 0; i < points; i++) {
		float length = (max_radius * 0.5) + (rand.Double() * max_radius * 0.5);

		float angle = (2 * M_PI / points) * i;

		float dx = cos(angle) * length;
		float dy = sin(angle) * length;
		nodes.push_back({ cx + dx, cy + dy });
	}

	closed_spline(nodes, parts, mindist, rand);

	// clamp and draw points
	for (auto &pnt: nodes) {
		pnt.x = std::min((double)rows - 1, pnt.x);
		pnt.x = std::max(0.0, pnt.x);
		pnt.y = std::min((double)rows - 1, pnt.y);
		pnt.y = std::max(0.0, pnt.y);
		put_point(bitset, pnt.x, pnt.y, rows, cols);
	}

	// close contour
	for(unsigned i = 0; i < nodes.size(); ++i) {
		line_between(bitset, nodes[i], nodes[(i + 1) % nodes.size()], rows, cols);
	}

	fill(bitset, rows, cols, cx, cy);
	create_patches(bitset, rows, cols, 15, 1, rand);
}
