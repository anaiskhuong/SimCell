#include "main/util.hpp"
#include "main/geom.hpp"
#include "main/cell.hpp"
#include "main/data.hpp"

extern int	*domin;
extern int	print, cycles, delays, hold, start;
extern FILE	*pdb, *msd, *rdf, *dom, *coord;
extern int	lastlook;

float accessArea ( Vec*, float*, int );
Vec shrinkWrap ( Vec*, int );
void rdfcalc ();
void domain ( int, float* );
int  setDomID ( Vec*, int*, int );


void looker ()
{
Cell	*world = Cell::world;
Data	*model = Data::model+world->model;
int	run = Data::frame;
int	n = world->kids, nn = n+1;
int	look = 0, waitime = 50;
int linkedCells = 0;

float	split;
	if (run <= delays) return; // files won't be open yet
	if (lastlook < waitime) return;

        if (run <= delays+hold+cycles+hold)
        {
            if (run%hold != 0) return;
//            Pi(run) Pi(lastlook) Pi(print) Pi(delays) NL
            lastlook = run;
            return;
            DO(i,n) { Cell *ci = world->child[i];
                        DO (j,ci->kids){ Cell *cj = ci->child[j];
                                if (cj->nlinks>0) {linkedCells = 1;}
                        }
//				fprintf(coord, "%d %d %d %f %f %f\n", run, ci->uid, linkedCells, ci->xyz.x, ci->xyz.y, ci->xyz.z);
//				fprintf(coord, "%f %f %f %d\n", ci->xyz.x, ci->xyz.y, ci->xyz.z, ci->uid);
                        linkedCells = 0;
                }
        }

	if (print > 0 ) // calc MSD and RDF plots
	{ float	d, dd, ddsum, ddRsum, ddGsum, ddXsum;
	  int	nR=0, nG=0, nX=0;
		Pt(Printing) Pi(run) NL
		ddsum = ddRsum = ddGsum = ddXsum = 0.0;
		// calc MSD over, all, Red, Green, eXcited (model==0 = G)
                DO(i,n) { Cell *ci = world->child[i];
                        dd = ci->vdata[0] || ci->xyz; // distance**2 from starting position
			ddsum += dd;
			if (ci->model==0) { ddGsum += dd; nG++; } else { ddRsum += dd; nR++; } // R/G
			if (ci->fdata[0]+ci->fdata[1] > 2.0+NOISE) { ddXsum += dd; nX++; } // excited
                }
		ddsum /= (float)n;
		if (nR) ddRsum /= (float)nR;
		if (nG) ddGsum /= (float)nG;
		if (nX) ddXsum /= (float)nX;
		Pr(ddsum) Pr(ddRsum) Pr(ddGsum) Pr(ddXsum) NL
		fprintf(msd,"%d   %f %f   %f %f   %f %f   %f %f\n", run-start,
			ddsum,sqrt(ddsum), ddRsum,sqrt(ddRsum), ddGsum,sqrt(ddGsum), ddXsum,sqrt(ddXsum));
		// excited cells
                DO(i,n) { Cell *ci = world->child[i];
			if (ci->idata[2] < 1) continue;
                        dd = ci->vdata[1] || ci->xyz; // distance**2 from starting position
			Pi(i) Pi(run) Pi(ci->idata[2]) Pr(dd) NL
		}
		// calc RDF
		rdfcalc();
	}
	if (run%200 == 0) // define clusters
	{ Vec	*doms = new Vec[nn],
		*bead = new Vec[nn*8];
	  float *temp = new float[nn];
	  int	*domid = new int[nn],
		*rank = new int[nn],
		maxin, max, id, in, im;
	  float pmax, pmin;
	  int	ndom;
		// define clusters
		domid[0] = 9999;
		// call dom() seeded by cell rank in each dimension (123=XYZ)
		domain(1,temp); DO1(i,n) doms[i].x = temp[i];
		domain(2,temp); DO1(i,n) doms[i].y = temp[i];
		domain(3,temp); DO1(i,n) doms[i].z = temp[i]*0.5;
		ndom = setDomID(doms,domid,n) - 1;
		DO1(i,n) { Cell *ci = world->child[i-1];
			ci->resn = domid[i]; // save value
		}
		printf("%3d clusters", ndom);
		sort(domid,rank,-nn); // -ve n = reverse sort
		maxin = 0;
		id = in = 0;
		DO(i,n) { int ri = rank[i], rid = domid[ri];
			if (rid == id) {
				if (look) printf("%3d ", ri);
				in++;
				if (in > maxin) { maxin = in; max = id; }
			} else { int type = world->child[ri-1]->model;
				if (in) domin[in]++;
				id++;
				if (look) printf("\n%3d %1d  =  %3d ", id,type,ri);
				in = 1;
			}
		} NL
		split = 0.0;
		DO1(k,ndom) { Vec out; float path, area, surf; char col;
			in = 0;
			DO1(i,n) { float domi = (float)domid[i]; Cell *ci = world->child[i-1]; char c;
				if (domid[i] != k) continue;
				if (ci->model) col = 'R'; else col = 'G';
				bead[in++] = ci->xyz; // include cell centre
				DO(j,8) bead[in++] = ci->child[j]->xyz;
			}
			out = shrinkWrap(bead,in);
			path = out.x; area = out.y; surf = out.z;
			if (k<3) split += path;
			in /= 9;
			if (look) printf("domain %d = %d %c cells: path = %f area = %f surf = %f dens = %f\n",
				k, in, col, path, area, surf, 10.0*(float)in/area);
			//fprintf(dom,"domain %d = %d %c cells: path = %f area = %f surf = %f dens = %f\n",
			//	k, in, col, path, area, surf, 10.0*(float)in/area);
		}
		Pi(run) Pr(split) NL
	}
	if (print > 1) { // dump coordinates (with domain id as res.No)
		print = 1;
		DO(i,n) { Cell *ci = world->child[i]; char c;
		  	if (ci->model) c = 'C'; else c = 'L';
			fprintf(pdb,"ATOM%7d  CA  %cYS A%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n", i,c,ci->resn,
				-ci->xyz.x,-ci->xyz.y,ci->xyz.z, 0.1*(float)ci->resn+10.0+10.0*ci->model);
		}
		fprintf(pdb,"TER\n");
		DO(i,n) { Cell *ci = world->child[i]; char c;
		  	if (ci->model) c = 'C'; else c = 'L';
			fprintf(pdb,"ATOM%7d  CA  %cYS B%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n", i,c,i+1,
				-ci->xyz.x,-ci->xyz.y,ci->xyz.z, 0.1*(float)ci->resn+10.0+10.0*ci->model);
		}
		fprintf(pdb,"TER\n");
		fprintf(pdb,"ATOM      0  CA  TRP C   0       0.000   0.000   0.000  0.00  0.00\n"); // set 0
	}
	if (print == -1) {
		Pt(Printing stopped) NL
		fclose(msd);
		fclose(rdf);
		fclose(pdb);
		fclose(dom);
		DO1(i,n) domin[i] = 0;
		print = 0;
	}
	if (print) lastlook = 0;
	sleep(1);
}

float accessArea ( Vec *a, float *seg, int n ) {
float	**d, *s, *t, perim, rad = 2.0, rrad = rad+rad;
int	*tag, incell = 9;
	tag = new int[n];
	s = new float[n];
	t = new float[n];
	d = new float*[n]; DO(i,n) d[i] = new float[n];
	DO(i,n) DO(j,n) d[i][j] =  a[i]|a[j];
	perim = 0.0;
	DO(i,n) { float sero, smin, tmax, ends;
		seg[i] = 0.0;
		DO(j,n) s[j] = t[j] = -99.9;
		DO(j,n)
		{ Vec v = (a[j]-a[i]).norm();
	  	  float alf, bet;
			if (i==j) continue;
			if (d[i][j] > rrad) continue;
			bet = angle2pi(v.x,v.y);	// bet = angle from +Y
			alf = acos(0.5*d[i][j]/rad);	// alf = angle to intersect
			s[j] = bet-alf;		// segment starting angle
			t[j] = bet+alf;		// segment terminal angle
		}
		// sort on start position
		sort(s,tag,-n);
		smin = 999.9;
		sero = tmax = -999.9;
		DO(j,n) { int k = tag[j];
			if (s[k] < -99.0) continue;
			if (sero < -999.0) {		// first contact
				sero = s[k];
			} else {			// check -ve first not lapped
				if (s[k] > sero+twoPI) break;
			}
			if (t[k] < tmax) continue;	// segment is fully covered
			if (s[k] < smin) smin = s[k];
			if (s[k] > tmax) {		// start a new segment
				smin = s[k];
				if (tmax > 0.0) {	// not the start
					seg[i] += smin-tmax;
				}
				tmax = -999.9;
			}
			if (t[k] > tmax) tmax = t[k];
			if (sero > 0.0) {	// check if past the start
				if (tmax>twoPI && tmax-twoPI > sero) break;
			} else {
				if (tmax > sero+twoPI) break;
			}
		}
		/*
		// 4 possible end configurations (but only one value)
		if (sero>0.0 && tmax<twoPI) ends = sero+twoPI-tmax; // sero + twoPI-tmax;
		if (sero>0.0 && tmax>twoPI) ends = sero+twoPI-tmax; // sero - (tmax-twoPI);
		if (sero<0.0 && tmax<twoPI) ends = sero+twoPI-tmax; // sero+twoPI - tmax;
		if (sero<0.0 && tmax>twoPI) // overlap
		*/
		ends = sero+twoPI-tmax;
		if (ends > 0.0) seg[i] += ends; // -ve ends = overlap
		perim += seg[i];
	}
	return perim;
}

float triArea ( Vec *a, int p, int q, int r ) {
	if (p<0 || q<0 || r<0) return 0.0;
	return ((a[q]-a[p])^(a[r]-a[p])).len()*0.5;
}

typedef struct {
	int *pole, n;
} List;

Vec shrinkWrap ( Vec *a, int n ) {
float	xmin,xmax, ymin,ymax;
int	start = 1, incell = 9;
float	**d, **w, area, perim, split,  cut = 2.0;
int	in, *path, *used;
List	*list = new List[4];
float	*surf, surface;
	DO(i,n) a[i].z = 0.0;
	DO(i,4) {
		list[i].n = 1;
		list[i].pole = new int[n];
	}
	used = new int[n]; DO(i,n) used[i] = 0;
	path = new int[n+1]; DO(i,n) path[i] = -1;
	surf = new float[n];
	surface = accessArea(a,surf,n);
	DO(i,n) if (surf[i] < NOISE) used[i] = -1;
	d = new float*[n]; DO(i,n) d[i] = new float[n];
	w = new float*[n]; DO(i,n) w[i] = new float[n];
	DO(i,n) DO(j,n)
	{ float x = a[i]|a[j]; // + randf()*0.1;
		d[i][j] = x;
		w[i][j] = exp(-(x*x)*0.1);
	}
	xmin = ymin =  9999.9;
	xmax = ymax = -9999.9;
	// find N.S.E.W poles
	DO(i,n) {
		if (i%incell==0) continue; // skip centroids
		if (a[i].y > ymax) { ymax = a[i].y; list[0].pole[0] = i; } 
		if (a[i].x > xmax) { xmax = a[i].x; list[1].pole[0] = i; } 
		if (a[i].y < ymin) { ymin = a[i].y; list[2].pole[0] = i; } 
		if (a[i].x < xmin) { xmin = a[i].x; list[3].pole[0] = i; } 
	}
	// eliminate any duplicate poles
	for (int i = 1; i<4; i++) if (list[0].pole[0]==list[i].pole[0]) list[i].n = -1;
	for (int i = 2; i<4; i++) if (list[1].pole[0]==list[i].pole[0]) list[i].n = -1;
	for (int i = 3; i<4; i++) if (list[2].pole[0]==list[i].pole[0]) list[i].n = -1;
	// find starting area
	area = 0.0;
	area += triArea(a,list[0].pole[0],list[1].pole[0],list[2].pole[0]);
	area += triArea(a,list[0].pole[0],list[3].pole[0],list[2].pole[0]);
	if (list[2].pole[0]<0) area = triArea(a,list[0].pole[0],list[1].pole[0],list[3].pole[0]);
	// find best points between poles
	LOOP
	{ int	got = 0,
		last = list[0].pole[0];
		for (int ii=3; ii>=0; ii--)
		{ List	*lii = list+ii;
		  int	len = lii->n;
			if (len < 0) continue;
			for (int i=len-1; i>=0; i--)
			{ float smax; int best, pole = lii->pole[i];
				used[pole] = used[last] = 1;
				// find best point j between last and pole
				in = 0;
				smax = 0.0;
				DO(j,n) { float score, ang, wR, wL;
					if (used[j]) continue;
					if (j==pole || j==last) continue;
					if (j%incell==0) continue; // skip centroids
					if (d[j][pole] > d[last][pole]) continue;
					if (d[j][last] > d[last][pole]) continue;
					ang = angle(a[last],a[j],a[pole]);
					if (ang < halfPI) continue;
					wR = wL = 1.0;
					DO(k,n) { Vec p,q,r; float e,f, wjk = w[j][k];
						if (start) wjk = 1.0; // consider all at the start
						if (d[j][k] > d[pole][last]) continue;
						if (pole==k || last==k || j==k) continue;
						p = a[pole]-a[j]; q = a[last]-a[j]; r = a[k]-a[j];
						e = (p.y*r.x - p.x*r.y)/(p.y*q.x - p.x*q.y);
						f = (q.y*r.x - q.x*r.y)/(q.y*p.x - q.x*p.y);
						if (k%incell==0) wjk *= 10.0; // up-weight centroids
						if (e>0.0 && f>0.0) wR += wjk; else wL += wjk;
					}
					score = ang + 10.0*surf[j];
					if (wR > wL) score /= wL; else score /= wR; // use smaller penalty
					if (score > smax) { int ok = 1;
						if (d[last][pole]<cut && ang<halfPI) ok = 0;
						if (ok) { best = j; smax = score; in++; }
					}
				}
				if (in) { float tri, sum = 0.0;
					DO(iii,4) {
						if (list[iii].n<0) continue;
						DO(jjj,list[iii].n)
						{ int	 apex = list[iii].pole[jjj];
						  float e,f, dpl = d[apex][pole]+d[apex][last];
						  Vec	p,q,r;
							if (apex==pole || apex==last) continue;
							p = a[pole]-a[apex]; q = a[last]-a[apex]; r = a[best]-a[apex];
							e = (p.y*r.x - p.x*r.y)/(p.y*q.x - p.x*q.y);
							f = (q.y*r.x - q.x*r.y)/(q.y*p.x - q.x*p.y);
							if (e>0.0 && f>0.0) { // clear view from apex
								 if (e+f<1.0) dpl = -dpl;
							} else { dpl = 0.0; }
							sum += dpl;
						}
					}
					tri = triArea(a,last,best,pole);
					if (sum > 0.0) area += tri; else area -= tri;
					got = 1;
					used[best] = 1;
					if (len>1) {	// insert
						in = 0;
						DO(k,lii->n-1)
						{ int	l, k0 = lii->pole[k], k1 = lii->pole[k+1];
							if (pole==k0 && last==k1) {
								in = 1;
								for (l=lii->n; l>k+1; l--) { // bump-up
									lii->pole[l] = lii->pole[l-1];
								}
								lii->pole[k+1] = best;
								break;
							}
						}
						if (in==0) lii->pole[len] = best;
					} else {	// add on end
						lii->pole[len] = best;
					}
					lii->n++;
				}
				last = pole;
			}
		}
		if (got==0) break;
		start = 0;
	}
	in = 0;
	DO(i,4) {
		if (list[i].n<0) continue;
		DO(j,list[i].n) path[in++] = list[i].pole[j];
	}
	split = perim = 0.0;
	DO(i,in) { float d,e,f;
		// if (path[i]>99) printf("%4d", path[i]); else printf("%3d", path[i]);
		if (i<1) continue;
		d = a[path[i-1]]|a[path[i]];
		perim += d;
		e = (a[path[i-1]]+a[path[i]]).len();	// double distance of segment from centre
		f = Data::model[0].sizes[0];		// diameter of the world
		if (e > f*0.6 ) continue;
		split += d;
	} // NL
/*
FILE *ut = fopen("test.pdb","w");
        DO(i,n) { char c;
		if (used[i] > 0) c = 'C'; else c = 'L';
		if (i%incell==0) c = 'Z';
                fprintf(ut,"ATOM%7d  CA  %cYS A%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n",
			i,c,i, -a[i].x, -a[i].y, a[i].z, surf[i]);
        }
fclose(ut);
	return Vec(perim,area,surface);
*/
	return Vec(split,area,surface);
}

void rdfcalc () {
// calculate the RDF for red/green cells
Cell	*world = Cell::world;
float	x,y,z;
float	allh, alls;
float	midx, midy, maxr, r;
float	smooth[999], base[999];
int	i,j, nred, ngrn, nall;
int	average = 0; // average over adjacent bins
int	count[999];
Vec	red[9999], grn[9999];
//
int	n=30;  // 20*20 = 400 - corners = 300 (100*pi)
float	scale = 33.3;
float	dense = 0.333; // 0.49999 max
float	c = scale*(float)(n/2);
float	xmin,xmax, ymin,ymax, zmin,zmax, size;
int	allbump, redbump, grnbump, rngbump;
int	bump = 4.50;
/*
long	seed;
	seed = (long)time(0);
	srand48(seed);
*/
	xmin = ymin = zmin = 9999.9;
	xmax = ymax = zmax = -9999.9;
	maxr = 0.5*world->len;
	// read in cells
	nred = ngrn = 0;
        DO(i,world->kids) { Cell *ci = world->child[i];
		if (ci->model==0) { grn[ngrn] = ci->xyz; ngrn++; }
		if (ci->model==1) { red[nred] = ci->xyz; nred++; }
		x = ci->xyz.x; y = ci->xyz.y; z = ci->xyz.z;
		if (x>xmax) xmax = x; if (x<xmin) xmin = x;
		if (y>ymax) ymax = y; if (y<ymin) ymin = y;
		if (z>zmax) zmax = z; if (z<zmin) zmin = z;
        }
	xmax = xmax-xmin;
	ymax = ymax-ymin;
	zmax = zmax-zmin;
	size = xmax;
	if (ymax>size) size = ymax;
	if (zmax>size) size = zmax;
        dense = PI*maxr*maxr*0.25/(float)(nred+ngrn);
        Pi(nred) Pi(ngrn) Pi(nred+ngrn) Pr(size) Pr(maxr) Pr(dense) NL
	for (i=0; i<100; i++) { float x = 0.5+(float)i, R = 0.5*size; //maxr;
		base[i] = (2.0*x/(PI*R*R))*(2.0*acos(x/(2.0*R))-(x/R)*sqrt(1.0-x*x/(4.0*R*R)));
	}
	// count all pairs rr, gg, rg
	allbump = redbump = grnbump = rngbump = 0;
	size = 50.0/maxr;
	size = 1.0;
	n = 0;
	for (i=0; i<999; i++) count[i] = 0;
	for (i=0; i<nred-1; i++) {
		for (j=i+1; j<nred; j++)
		{ float d = vdif(red[i],red[j]);
		  int	bin = (int)(d*size);
			count[bin]++; n++;
			if (d<bump) allbump++;
		}
	}
	for (i=0; i<ngrn-1; i++) {
		for (j=i+1; j<ngrn; j++)
		{ float d = vdif(grn[i],grn[j]);
		  int	bin = (int)(d*size);
			count[bin]++; n++;
			if (d<bump) allbump++;
		}
	}
	for (i=0; i<nred; i++) {
		for (j=0; j<ngrn; j++)
		{ float d = vdif(red[i],grn[j]);
		  int	bin = (int)(d*size);
			count[bin]++; n++;
			if (d<bump) allbump++;
		}
	}
	for (i=0; i<100; i++)
	{ float	counti = (float)count[i], fx, x=0.5+(float)i;
		if (average) {
			if (i==0 || i==99) { counti = 0.0; continue; }
			counti = 0.25*(float)(count[i-1]+count[i]+count[i]+count[i+1]);
		}
		fx = (float)counti/(float)n;
		fprintf(rdf,"ALL %f %f %f %f\n",x,fx,fx/base[i],base[i]);
	}
	fprintf(rdf,"ALL\n");
	// count red pairs
	n = 0;
	for (i=0; i<999; i++) count[i] = 0;
	for (i=0; i<nred-1; i++) {
		for (j=i+1; j<nred; j++)
		{ float d = vdif(red[i],red[j]);
		  int	bin = (int)(d*size);
			count[bin]++; n++;
			if (d<bump) redbump++;
		}
	}
	for (i=0; i<100; i++)
	{ float	counti = (float)count[i], fx, x=0.5+(float)i;
		fx = (float)counti/(float)n;
		fprintf(rdf,"RED %f %f %f\n",x,fx,fx/base[i]);
	}
	fprintf(rdf,"RED\n");
	// count green pairs
	n = 0;
	for (i=0; i<999; i++) count[i] = 0;
	for (i=0; i<ngrn-1; i++) {
		for (j=i+1; j<ngrn; j++)
		{ float d = vdif(grn[i],grn[j]);
		  int	bin = (int)(d*size);
			count[bin]++; n++;
			if (d<bump) grnbump++;
		}
	}
	for (i=0; i<100; i++)
	{ float	counti = (float)count[i], fx, x=0.5+(float)i;
		fx = (float)counti/(float)n;
		fprintf(rdf,"GRN %f %f %f\n",x,fx,fx/base[i]);
	}
	fprintf(rdf,"GRN\n");
	// count red/green pairs
	n = 0;
	for (i=0; i<999; i++) count[i] = 0;
	for (i=0; i<nred; i++) {
		for (j=0; j<ngrn; j++)
		{ float d = vdif(red[i],grn[j]);
		  int	bin = (int)(d*size);
			count[bin]++; n++;
			if (d<bump) rngbump++;
		}
	}
	for (i=0; i<100; i++)
	{ float	counti = (float)count[i], fx, x=0.5+(float)i;
		fx = (float)counti/(float)n;
		fprintf(rdf,"RNG %f %f %f\n",x,fx,fx/base[i]);
	}
	fprintf(rdf,"RNG\n");
	Pi(allbump) Pi(redbump) Pi(grnbump) Pi(rngbump) NL
}
