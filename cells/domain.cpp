#include "main/util.hpp"
#include "main/geom.hpp"
#include "main/cell.hpp"
#include "main/data.hpp"

#define MINDOM 1
#define MAXDOM 1000

typedef struct {
	char	*res;
	float	*acc;
	int	*dom;
	int	*rid;
	Vec	*ca;
	int	len;
} Seq;

int  pairup ( float**, float*, float, int );
void setSij ( float**, float*, float, int );
void evolve ( float**, float*, int );
int  setDomID ( Vec*, int*, int );

void domain ( int dim, float *doms )
{
Cell	*world = Cell::world;
int	n = world->kids, nn = n+1; 
float	**mat = new float*[nn]; DO1(i,n) mat[i] = new float[nn];
Seq	*seq = new Seq;
float	size = 4.0; // distance at which cell are in contact
	seq->res = new char[nn];
	seq->acc = new float[nn];
	seq->dom = new int[nn];
	seq->rid = new int[nn];
	seq->ca  = new Vec[nn];
	// copy cell data to domain structure
	seq->len = n;
	DO1(i,n) { Cell *ci = world->child[i-1];
		if (ci->model) {
			seq->res[i] = 'C';
			seq->acc[i] =  1.0;
		} else {
			seq->res[i] = 'L';
			seq->acc[i] = -1.0;
		}
		seq->rid[i] = i;
		if (dim) seq->dom[i] = ci->ranks[dim-1];   // set id to rank in XYZ
		    else seq->dom[i] = (int)(doms[i]+0.5); // set id to consensus for final pass
		if (ci->model) {
			seq->dom[i] = -seq->dom[i] - 10;
		} else {	// separate cell types (adding +/-10 to deal with ambiguous rank 0) 
			seq->dom[i] =  seq->dom[i] + 10;
		}
		seq->acc[i] = (float)seq->dom[i]; // acc[] initialises dom[] labels
		seq->ca[i] = ci->xyz;
	}
	// define domains
	DO1(i,n) DO1(j,n) mat[i][j] = seq->ca[i] | seq->ca[j];
	pairup(mat,seq->acc,size,n);	// equalise values for cell in contact
	setSij(mat,seq->acc,size,n);	// convert mat to coupling strength (with screening)
	evolve(mat,seq->acc,n);		// iterate to consensus domain values
	DO1(i,n) DO1(j,n) mat[i][j] = seq->ca[i] | seq->ca[j];
	pairup(mat,seq->acc,size,n);	// re-equalise contact pair values
	DO1(i,n) doms[i] = seq->acc[i];
}

int pairup ( float **mat, float *dom, float size, int n )
{ // equalise value for cells in contact (closer than <diamet>er)
float	*mod = new float [n+1];
int	changes;
	DO1(i,n) mod[i] = dom[i];
	DO1(k,n) { float sum;
		changes = 0;
		DO1(j,n) { // loop over cells
			sum = 0.0;
			DO1(i,n) { // loop over neighbours
				if (i==j) continue;
				if (dom[i]*dom[j]<0.0) continue;
				if (mat[i][j] > size) continue;
				// force to values of touching pair to be the closer to zero
				if (dom[j] > 0.0) {
					if (dom[i] < dom[j]) {
						mod[j] = dom[i]; // close neighbour i is lower 
						changes++;
					}
				} else {
					if (dom[i] > dom[j]) {
						mod[j] = dom[i]; // close neighbour i is higher
						changes++;
					}
				}
			}
		}
		if (changes==0) break;
		DO1(j,n) dom[j] = mod[j];
	}
}

void setSij ( float **mat, float *ave, float size, int n )
{
int	i, j, k, nn = n+1;
float	spread = 7.0;
float	rad = 0.5*size,
	r = rad*1.0; // lower rad factor -> faster drop in screening with increasing h
	DO1(i,n) { float s, hh, x;
		DO1(j,n) { float a, b, c = mat[i][j];
			if (i>=j) continue;
			if (c > spread) continue;
			if (ave[i]*ave[j]<0.0) continue; // -ve = diff
			// check like-cells for interveaning unlike cell
			DO1(k,n) {
				if (k<i) a = mat[k][i]; else a = mat[i][k];
				if (k<j) b = mat[k][j]; else b = mat[j][k];
				if (ave[i]*ave[k]>0.0) continue;
				if (a > c || b > c) continue;
				x = (b*b-a*a+c*c)/(2.0*c);
				hh = b*b - x*x; 
				s = 1.0+exp(-hh/(r*r));
				mat[j][i] *= s; // modify distance in [j][i] half
			}
		}
	}
	DO1(i,n) DO1(j,n) {
		if (i>=j) continue;
		mat[i][j] = mat[j][i]; // copy back changes to [i][j] half
	}
	DO1(i,n) {
		DO1(j,n) { float d = mat[i][j];
			if (i==j) continue;
			mat[i][j] = 1.0-1.0/(1+exp((5.0-d)*2.0)); // switch function (mean 5, slope 3)
			if (d>spread) mat[i][j] = -1.0; // -ve = out of range
				else  mat[i][i] += 1.0; // count in-range on diag.
			if (ave[i]*ave[j]<0.0) mat[i][j] = -1.0; // skip mixed types
		}
	}
}

void evolve ( float **mat, float *ave, int n )
{
int	i, j, k;
float	*tmp, rms, sum, shift, step, mean, wt;
int	in, jumps, **list;
float	*was, *now;
int	cycles = n*5;
	was = new float[n+1];
	now = new float[n+1];
	for (j=1; j<=n; j++) was[j] = now[j] = ave[j];
	list = new int*[n+1];
	for (i=1; i<=n; i++) { // mat[diag] used to count neighbours
		list[i] = new int[(int)(mat[i][i]+1.5)];
		k = 1;
		for (j=1; j<=n; j++) { // fill <list[]> with neighbours (mat>0)
			if (i==j) continue;
			if (mat[i][j] < 0.0) continue;
			list[i][k] = j;
			k++;
		}
		list[i][0] = k; // hold <list> length in 0
	}
	// shift values towards neighbours
	shift = 1.0;
	step = 2.0/(float)cycles;
	DO1(k,cycles) { int ii;
		DO1(j,n) { // loop over cells
			sum = 0.0;
			for (ii=1; ii<list[j][0]; ii++) { // loop over neighbours
				i = list[j][ii];
				// gather consensus for shift
				if (was[i]>was[j]) sum += mat[i][j]; // vote for bigger value on j
				if (was[i]<was[j]) sum -= mat[i][j]; // vote for smaller value
			}
			if (sum>0.0) now[j] += shift;
			if (sum<0.0) now[j] -= shift;
		}
		rms = 0.0;
		DO1(j,n) { float a;
			a = ave[j];
			ave[j] = 0.5*(now[j]+was[j]);
			a -= ave[j];
			rms += a*a;
			was[j] = now[j];
		}
		rms = sqrt(rms/(float)n);
		if (rms<NOISE || k>cycles/2) shift -= step; 
		if (shift < 0.0) break;
	}
}

int setDomID ( Vec *doms, int *domid, int n ) {
// assign integer ID values based on XYZ seeded float values
// NB struct { int a,b; float s; char c; } Pairs (in geom.hpp)
Cell	*world = Cell::world;
Pairs	*pair = new Pairs [(n*n-n)/2];
int	*rank = new int [(n*n-n)/2];
float	d, cut = 1.0;
int	id, m = 0;
	DO1(i,n) { Cell *ci = world->child[i-1];
		domid[i] = 0;
		DO1(j,n) { Cell *cj = world->child[j-1];
			if (j >= i ) continue;
			if (ci->model != cj->model) continue;
			if (doms[i].x>0 && doms[j].y<0)
				{ Pt(XYerror) Pi(i) Pi(j) Pr(doms[i].x) Pr(doms[j].y) NL exit(1); }
			if (doms[i].x<0 && doms[j].y>0)
				{ Pt(YXerror) Pi(i) Pi(j) Pr(doms[i].x) Pr(doms[j].y) NL exit(1); }
			d = doms[i] | doms[j]; // distance between XYZ label values
			if (d > cut) continue;
			pair[m].a = i; pair[m].b = j; pair[m].s = d;
			if (ci->model) pair[m].c = 'R'; else pair[m].c = 'G';
			m++;
		}
	}
	sort(pair,rank,-m); // -m makes reverse sort
	/*
	DO(i,m) { int ri = rank[i];
		Pi(i) Pi(ri) Pi(pair[ri].a) Pi(pair[ri].b) Pc(pair[ri].c) Pr(pair[ri].s) NL
	}
	*/
	id = 1;
	DO(i,m) { int	ri = rank[i], ai = pair[ri].a, bi = pair[ri].b;
		if (domid[ai]==0 && domid[bi]==0) {
			domid[ai] = domid[bi] = id; // set closest free pair to next id
		} else { continue; }
		DO1(j,n) { int got = 0;
			DO(k,m) { int	rk = rank[k], ak = pair[rk].a, bk = pair[rk].b;
				if (domid[ak] && domid[bk]) continue;
				if (domid[ak]) { domid[bk] = id; got = 1; break; }
				if (domid[bk]) { domid[ak] = id; got = 1; break; }
			}
			if (got==0) break; // no more found
		}
		id++;
	}
	DO1(i,n) { // assign Jonny-no-mates
		if (domid[i]==0) domid[i] = id++;
	}
	/*
	DO1(i,n) {
		Pt(DOMS)
		printf("%4d %4d %9.3f %9.3f %9.3f\n", i, domid[i], doms[i].x, doms[i].y, doms[i].z);
	//	doms[i].x+randf(), doms[i].y+randf(), doms[i].z+randf());
	} NL
	*/
	return id;
			
}
