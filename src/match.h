#ifndef _MINIMUM_WEIGHT_PERFECT_MATCHING_H_
#define _MINIMUM_WEIGHT_PERFECT_MATCHING_H_

void MinimumWeightPerfectMatching(
	// number of vertices
	int N,
	// Returns cost for edge i-j
	double (*cost)(int i, int j, void* data),
	// Returned matching
	int *matches, // length N, so each pair is represented twice
	// matches[0] gives the vertex index matched with vertex 0,
	// therefore matches[matches[i]] == i should always be true
	void *data
);

#endif // _MINIMUM_WEIGHT_PERFECT_MATCHING_H_
