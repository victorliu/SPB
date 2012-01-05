#include <memory>
#include <cstring>
#include <float.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>

void match(
	int N,
	double *total_weight,
	int *nmatch,
	double (*cost)(int i, int j, void* data), void* data,
	int *basis,
	int *mem,
	int *ka,
	int *kb,
	int *sm,
	int *tma,
	int *tmb,
	int *m1,
	double *y1,
	double *y2,
	double *dplus,
	double *dminus);

void scan1(
	int nb1,
	int N,
	double (*cost)(int i, int j, void* data), void* data,
	int *basis,
	int *mem,
	int *ka,
	int *kb,
	int *sm,
	int *tma,
	int *tmb,
	double *y1,
	double *y2,
	double *dplus,
	double *dminus,
	int *m1);

void scan2(
	int nb,
	int N,
	double (*cost)(int i, int j, void* data), void* data,
	int *basis,
	int *mem,
	int *ka,
	int *kb,
	int *sm,
	int *tma,
	int *tmb,
	double *y1,
	double *y2,
	double *dplus,
	double *dminus);

void MinimumWeightPerfectMatching(
	int N,
	double (*cost)(int i, int j, void* data),
	int *matches, // length N
	void* data
){
	double *dtemp = (double*)malloc(sizeof(double)*4*N);
	double *y1 = dtemp;
	double *y2 = y1 + N;
	double *dplus = y2 + N;
	double *dminus = dplus + N;
	int *itemp = (int*)malloc(sizeof(int)*8*N);
	int *basis = itemp;
	int *mem = basis + N;
	int *ka = mem + N;
	int *kb = ka + N;
	int *sm = kb + N;
	int *tma = sm + N;
	int *tmb = tma + N;
	int *m1 = tmb + N;
	memset(dtemp, 0, sizeof(double)*4*N);
	memset(itemp, 0, sizeof(int)*8*N);
	double total_weight;
	match(N, &total_weight, matches, cost, data, basis, mem, ka, kb, sm, tma, tmb, m1, y1, y2, dplus, dminus);
	free(dtemp);
	free(itemp);
}


void match(
	int N, // number of nodes (must be even)
	double *total_weight, // cost of optimal matching
	int *nmatch, // matching array
	double (*cost)(int i, int j, void* data), void *data,
	// The following are work arrays each of length N
	int *basis,
	int *mem,
	int *ka,
	int *kb,
	int *sm,
	int *tma,
	int *tmb,
	int *m1,
	double *y1,
	double *y2,
	double *dplus,
	double *dminus
){
	const int top = N + 1;
	for(int n1 = 0; n1 < N; ++n1){
		basis[n1] = n1;
		mem[n1] = n1;
		y1[n1] = 0.0;
		y2[n1] = 0.0;
		sm[n1] = top;
		tma[n1] = top;
		tmb[n1] = top;
		nmatch[n1] = top;
		dplus[n1] = DBL_MAX;
		dminus[n1] = DBL_MAX;
		ka[n1] = -1;
		kb[n1] = n1;
	}

	// Start of main procedure
	for(int n1 = 0; n1 < N; ++n1){
		if(nmatch[n1] != top){ continue; }
		int nn = -1;
		double d = DBL_MAX;
		for(int n2 = 0; n2 < N; ++n2) {
			if(n1 == n2){ continue; }
			double newcost = cost(n1,n2,data) - y1[n2];
			if(newcost < d){
				d = newcost;
				nn = -1;
				if(nmatch[n2] == top){
					nn = n2;
				}
			}else if(newcost == d && nn < 0 && nmatch[n2] == top){
				nn = n2;
			}
		}
		if(nn >= 0){
			y1[n1] = d;
			nmatch[n1] = nn;
			nmatch[nn] = n1;
		}
	}
	

	// Initial labeling
	int nn = 0;
	for(int ni = 0; ni < N; ++ni){
		if(nmatch[ni] != top){ continue; }
		nn++;
		sm[ni] = -1;
		dplus[ni] = 0.0;
		double Y1B = y1[ni];
		for(int nk = 0; nk < N; ++nk){
			if(ni == nk){ continue; }
			double newcost = cost(ni,nk,data) - Y1B - y1[nk];
			if(newcost < dminus[nk]){
				dminus[nk] = newcost;
				ka[nk] = ni;
			}
		}
	}
	if(nn <= 1){ goto finalize; }

	// Examination of the labeling and decision for the next step
	{
make_decision:
		int n1, n2;
		int nka, nkb;
		double dbest = DBL_MAX;
		int nbest;
		for(int nb = 0; nb < N; ++nb) {
			if(basis[nb] != nb){ continue; }
			double d = dminus[nb];
			if(sm[nb] >=  top){
				if(tma[nb] >= top){
					if(d >= dbest){ continue; }
					nbest = nb;
					dbest = d;
				}
				if(mem[nb] != nb){
					d += y1[nb];
					if(d < dbest){
						nbest = nb;
						dbest = d;
					}
				}
			}else{
				d = 0.5*(d + dplus[nb]);
				if(d <= dbest){
					nbest = nb;
					dbest = d;
				}
			}
		}

		if(tma[nbest] < top){
			// Expansion of a T labeled blossom
			n1 = mem[nbest];
			int nb3 = n1;
			nka = ka[n1];
			{
				int nk2 = n1;
				do{
					int nk1 = nk2;
					nkb = kb[nk1];
					double y1b = y1[nk1];

					do{
						basis[nk2] = nk1;
						y2[nk2] -= y1b;
						if(nk2 == nkb){ break; }
						nk2 = mem[nk2];
					}while(1);
					nk2 = mem[nkb];
					mem[nkb] = nk1;
				}while(nk2 != nka);
			}
			{
				double y1b = dplus[n1];
				y1[nbest] = y1b;
				mem[nbest] = nka;
				int nk2 = nka;
				do{
					y2[nk2] -= y1b;
					if(nk2 == nbest){ break; }
					nk2 = mem[nk2];
				}while(1);
			}

			int nk1 = nmatch[nbest];
			int nb = basis[sm[basis[nk1]]];
			if(nb != nbest){
				{
					int nb2 = nb;
					int nk;
					do{
						nk = tma[nb2];
						int nb1 = basis[nk];
						if(nb1 == nbest){break; }
						nb2 = sm[nb1];
						nb2 = basis[nb2];
					}while(1);

					tma[nb] = tma[nbest];
					tma[nbest] = tmb[nb2];
					tmb[nb] = tmb[nbest];
					tmb[nbest] = nk;
				}
				int nk3 = sm[nb];
				nb3 = basis[nk3];
				int nk4 = sm[nb3];
				sm[nb] = top;
				nmatch[nb] = nk1;
				int nb1 = nb3;
				do{
					nk1 = tma[nb1];
					int nk2 = tmb[nb1];
					tma[nb1] = nk4;
					tmb[nb1] = nk3;
					sm[nb1] = nk1;
					nmatch[nb1] = nk1;
					int nb2 = basis[nk1];
					nmatch[nb2] = nk2;
					nk3 = sm[nb2];
					sm[nb2] = nk2;
					if(nb2 == nbest){ break; }
					nb1 = basis[nk3];
					nk4 = sm[nb1];
					tma[nb2] = nk3;
					tmb[nb2] = nk4;
				}while(1);
			}
			int nk2 = tmb[nb];
			int nb1 = basis[nk2];
			dminus[nb1] = dbest;
			n1 = -1;
			if(nb1 != nb){
				nk1 = tma[nb1];
				nb3 = basis[nk1];
				tma[nb1] = tma[nb];
				tmb[nb1] = nk2;

				int nk;
				do{
					nk = sm[nb1];
					sm[nb1] = top;
					int nb2 = basis[nk];
					nk = tma[nb2];
					tma[nb2] = top;
					n2 = tmb[nb2];
					tmb[nb2] = n1;
					n1 = nb2;
					dplus[nb2] = dbest;
					nb1 = basis[nk];
					dminus[nb1] = dbest;
				}while(nb1 != nb);
				tma[nb] = n2;
				tmb[nb] = nk;
				sm[nb] = top;

				if(nb3 == nb){ goto do_scan1; }
			}
			nb1 = -1;
			{
				int nb2 = nb3;
				do{
					int nk = sm[nb2];
					sm[nb2] = top;
					tma[nb2] = top;
					tmb[nb2] = nb1;
					nb1 = basis[nk];
					
					nk = tma[nb1];
					sm[nb1] = top;
					tma[nb1] = top;
					tmb[nb1] = nb2;
					nb2 = basis[nk];
				}while(nb2 != nb);
			}

			scan2(nb1, N, cost, data, basis, mem, ka, kb, sm, tma, tmb,
				y1, y2, dplus, dminus);
do_scan1:
			while(n1 >= 0){
				int nb = n1;

				scan1(nb, N, cost, data, basis, mem, ka, kb, sm, tma, tmb,
					y1, y2, dplus, dminus, m1);

				n1 = tmb[nb];
				tmb[nb] = top;
			}
		}else if(sm[nbest] >= top){
			// Growing an alternating tree by adding two edges

			tma[nbest] = ka[nbest];
			tmb[nbest] = kb[nbest];
			int nmb = basis[nmatch[nbest]];
			dplus[nmb] = dbest;
			sm[nmb] = nmatch[nmb];

			scan1(nmb, N, cost, data, basis, mem, ka, kb, sm, tma, tmb,
				y1, y2, dplus, dminus, m1);
		}else{
			nka = ka[nbest];
			nkb = kb[nbest];
			n1 = nbest;
			int nb1 = n1;
			n2 = basis[nka];
			int nb2 = n2;
			do{
				tma[nb1] = nb2;
				int nk = sm[nb1];
				if(nk < 0){ break; }
				nb2 = basis[nk];
				nb1 = tma[nb2];
				nb1 = basis[nb1];
			}while(1);
			int nb = nb1;
			nb1 = n2;
			nb2 = n1;
			while(tma[nb1] >= top){
				tma[nb1] = nb2;
				int nk = sm[nb1];
				if(nk < 0){ goto augment_matching; }
				nb2 = basis[nk];
				nb1 = tma[nb2];
				nb1 = basis[nb1];
			}
			while(nb1 != nb){
				int nk = tma[nb];
				tma[nb] = top;
				nb = basis[nmatch[nk]];
			}
			
			// Shinking a blossom
			double yb = y1[nb] + dbest - dplus[nb];
			y1[nb] = 0.0;
			int nk1 = nb;
			do{
				y2[nk1] += yb;
				nk1 = mem[nk1];
			}while(nk1 != nb);
			int memsave = mem[nb];
			if(nb == n2){
				n2 = n1;
				nb2 = tma[nb];
			}
			do{
				mem[nk1] = nb2;
				int nm = nmatch[nb2];
				sm[nb2] = nm;
				double y1b = y1[nb2] + dminus[nb2] - dbest;
				nk1 = nb2;
				int nk2;
				do{
					nk2 = nk1;
					y2[nk2] += y1b;
					basis[nk2] = nb;
					nk1 = mem[nk2];
				}while(nk1 != nb2);
				kb[nb2] = nk2;
				y1[nb2] = y1b;
				nb1 = basis[nm];
				mem[nk2] = nb1;
				y1b = y1[nb1] + dbest - dplus[nb1];
				nk2 = nb1;

				do{
					nk1 = nk2;
					y2[nk1] += y1b;
					basis[nk1] = nb;
					nk2 = mem[nk1];
				}while(nk2 != nb1);
				kb[nb1] = nk1;
				y1[nb1] = y1b;
				if(n2 != nb1){
					nb2 = tma[nb1];
					tma[nb1] = tmb[nb2];
					tmb[nb1] = tma[nb2];
					continue;
				}
				
				if(n2 != nbest){
					tma[n2] = nkb;
					tmb[n2] = nka;
					if(nb != nbest){
						n2 = n1;
						nb2 = tma[nb];
						continue;
					}
					break;
				}else{
					tma[nbest] = nka;
					tmb[nbest] = nkb;
					break;
				}
			}while(1);

			mem[nk1] = memsave;
			n1 = mem[nb];
			ka[n1] = memsave;
			dplus[n1] = yb;
			tma[nb] = top;
			dplus[nb] = dbest;

			scan1(nb, N, cost, data, basis, mem, ka, kb, sm, tma, tmb,
				y1, y2, dplus, dminus, m1);
		}
		goto make_decision;

		// Augmentation of matching
		// Exchange of the matching and non matching edges along the augmenting path
augment_matching:
		{
			int nb = n1;
			int nk = nka;
			do{
				int nb1 = nb;
				do{
					nmatch[nb1] = nk;
					nk = sm[nb1];
					tma[nb1] = top;
					if(nk < 0){ break; }
					int nb2 = basis[nk];
					int nk1 = tma[nb2];
					nk = tmb[nb2];
					nb1 = basis[nk1];
					nmatch[nb2] = nk1;
				}while(1);
				if(nb != n1){ break; }
				nb = n2;
				nk = nkb;
			}while(1);
		}

		// Removing all labels on non exposed base nodes
		for(int nb = 0; nb < N; ++nb){
			if(basis[nb] != nb){ continue; }
			if(sm[nb] >= top){
				if(tma[nb] < top){
					double d = dminus[nb] - dbest;
					y1[nb] += d;
					tma[nb] = top;
					tmb[nb] = top;
				}
			}else{
				double d = dbest - dplus[nb];
				y1[nb] += d;
				sm[nb] = top;
				if(nmatch[nb] != top){
					dplus[nb] = DBL_MAX;
				}else{
					sm[nb] = -1;
					dplus[nb] = 0.0;
				}
			}
			dminus[nb] = DBL_MAX;
		}
		nn -= 2;
		if(nn <= 1){ goto finalize; }

		// Determination of the new dminus values
		for(int n1 = 0; n1 < N; ++n1){
			int nb1 = basis[n1];
			if(sm[nb1] >= 0){ continue; }
			double y1b = y1[nb1];
			double y2b = y2[n1];

			for(int n2 = 0; n2 < N; ++n2) {
				int nb2 = basis[n2];
				if(nb1 == nb2){ continue; }
				if(n1 == n2){ continue; }
				double newcost = cost(n1,n2,data) - y1b - y2b - y1[nb2] - y2[n2];
				if(newcost < dminus[nb2]){
					ka[nb2] = n1;
					kb[nb2] = n2;
					dminus[nb2] = newcost;
				}
			}
		}
		goto make_decision;
	}
	// Generation of the original graph by expansion of all shrunken blossoms
finalize:

	*total_weight = 0;
	for(int nb1 = 0; nb1 < N; ++nb1) {
		if(basis[nb1] != nb1 || sm[nb1] < -1){ continue; }
		int n2 = nmatch[nb1];
		int nb2 = basis[n2];
		int n1 = nmatch[nb2];
		sm[nb1] = -2;
		sm[nb2] = -2;
		if(n1 == n2){ continue; }
		double nc = cost(n1,n2,data);
		//double d = nc - y1[nb1] - y1[nb2] - y2[n1] - y2[n2];
		*total_weight += nc;
	}
	for(int n1 = 0; n1 < N; ++n1) {
restart_loop:
		int nb = basis[n1];
		if(nb == n1){ continue; }
		int nk2 = mem[nb];
		int nka = ka[nk2];
		int nb3 = nk2;
		double yb = dplus[nk2];
		do{
			int nk1 = nk2;
			int nkb = kb[nk1];
			double y1b = y1[nk1];
			do{
				basis[nk2] = nk1;
				y2[nk2] -= y1b;
				if(nk2 == nkb){ break;}
				nk2 = mem[nk2];
			}while(1);
			nk2 = mem[nkb];
			mem[nkb] = nk1;
		}while(nk2 != nka);
		y1[nb] = yb;
		mem[nb] = nka;
		nk2 = nka;
		do{
			y2[nk2] -= yb;
			if(nk2 == nb){ break; }
			nk2 = mem[nk2];
		}while(1);

		int nk = nmatch[nb];
		int nk1 = basis[nk];
		nk1 = nmatch[nk1];
		int nb1 = basis[nk1];
		if(nb == nb1){
			goto skip_forward;
		}
		nmatch[nb1] = nk;
		nb3 = tma[nb1];
		nb3 = basis[nb3];
		do{
			int nb2 = basis[sm[nb1]];
			int nk1 = tma[nb2];
			nk2 = tmb[nb2];
			nb1 = basis[nk1];
			nmatch[nb1] = nk2;
			nmatch[nb2] = nk1;
			if(nk1 == nk2){
				goto restart_loop;
			}
			double nc = cost(nk1,nk2,data);
			double d = nc - y1[nb1] - y1[nb2] - y2[nk1] - y2[nk2];

			if(fabs(d) > DBL_EPSILON){
				fprintf(stderr,
					"Optimality conditions are violated at edge %3d <--> %3d)\n",
					nk1, nk2);
			}

			*total_weight += nc;
		}while(nb1 != nb);
loop:
		if(nb3 == nb)  goto restart_loop;
skip_forward:
		int n2 = sm[nb3];
		int nb2 = basis[n2];
		int n3 = sm[nb2];
		if(n2 == n3){ goto restart_loop; }
		
		{
			double nc = cost(n2,n3,data);
			double d = nc - y1[nb2] - y1[nb3] - y2[n2] - y2[n3];

			if(fabs(d) > DBL_EPSILON){
				fprintf(stderr,
					"Optimality conditions are violated at edge %3d <--> %3d)\n",
					n2, n3);
			}

			*total_weight += nc;
		}
		n3 = tma[nb2];
		nb3 = basis[n3];
		goto loop;
	}
}

void scan1(
	int nb1,
	int N,
	double (*cost)(int i, int j, void *data), void *data,
	int *basis,
	int *mem,
	int *ka,
	int *kb,
	int *sm,
	int *tma,
	int *tmb,
	double *y1,
	double *y2,
	double *dplus,
	double *dminus,
	int *m1
){ // Scanning of node nb1
	const int top = N + 1;
	double d1 = dplus[nb1] - y1[nb1];
	dminus[nb1] = DBL_MAX;
	double d2 = d1 - y2[nb1];
	tma[nb1] = -1;
	int i1 = 0;
	for(int n2 = 0; n2 < N; ++n2) {
		int nb2 = basis[n2];
		if(tma[nb2] < top){ continue; }
		m1[i1] = n2;
		i1++;
		if(nb1 == n2){ continue; }
		double newcost = cost(nb1,n2,data) + d2 - y1[nb2] - y2[n2];
		if(newcost < dminus[nb2]){
			ka[nb2] = nb1;
			kb[nb2] = n2;
			dminus[nb2] = newcost;
		}
	}
	tma[nb1] = top;
	int n1 = mem[nb1];
	while(n1 != nb1){
		d2 = d1 - y2[n1];
		for(int i2 = 0; i2 < i1; ++i2) {
			int n2 = m1[i2];
			int nb2 = basis[n2];
			if(n1 == n2){ continue; }
			double newcost = cost(n1,n2,data) + d2 - y1[nb2] - y2[n2];
			if(newcost < dminus[nb2]){
				ka[nb2] = n1;
				kb[nb2] = n2;
				dminus[nb2] = newcost;
			}
		}
		n1 = mem[n1];
	}
	return;
}

void scan2(
	int nb,
	int N,
	double (*cost)(int i, int j, void *data), void *data,
	int *basis,
	int *mem,
	int *ka,
	int *kb,
	int *sm,
	int *tma,
	int *tmb,
	double *y1,
	double *y2,
	double *dplus,
	double *dminus
){
	// Scanning of node nb
	const int top = N + 1;
	do{
		int nb1 = nb;
		nb = tmb[nb1];
		tmb[nb1] = top;
		double d = DBL_MAX;
		int nka = -1;
		int nkb = -1;
		int n1 = nb1;
		double y1b = y1[nb1];
		do{
			double y2b = y2[n1];
			int n2;
			for(n2 = 0; n2 < N; ++n2) {
				int nb2 = basis[n2];
				if(sm[nb2] >= top){ continue; }
				if(n1 == n2){ continue; }
				double newcost = cost(n1,n2,data) - y1b - y2b - y1[nb2] - y2[n2] + dplus[nb2];
				if(newcost < d){
					nka = n2;
					nkb = n1;
					d = newcost;
				}
			}
			n1 = mem[n1];
		}while(n1 != nb1);
		ka[nb1] = nka;
		kb[nb1] = nkb;
		dminus[nb1] = d;
	}while(nb >= 0);
}
