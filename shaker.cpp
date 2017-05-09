#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

void shakeCell ( Cell*, Vec );
void spinsCell (Cell* );
void patchCell ( Cell*, float );
float fixgeom( Cell*, float );
float fixdist( Cell*, float );

void shaker () {
Vec	disp;
Cell	*world = Cell::world;
	DO(i,world->kids) // don't shake the world, just its contents
	{ Cell	*child = world->child[i];
		if (child->live < 1) continue;
		if (child->parent->solid > 0) continue;
		if (child->empty) continue;
		spinsCell(child);
		shakeCell(child,disp);
	}
}

void shakeCell ( Cell *cell, Vec disp )
{
float	kick;
int	i, m, n, id, level;
Data	*param = Data::model+cell->model;
int	moltype = param->moltype;
	if (cell->empty) return;
	level = cell->level;
	kick = param->kicks[level];
	if (kick < 0.0) kick = -kick; // -ve = no spin
//	if (moltype==1 && level==depth) { Cells *pa = cell->parent;
//		if ( pa->type==2 && pa->sort==1) kick *= 0.1; // keep calm in an RNA stem
//	}
	if (kick > NOISE) { // don't apply zero kick
		if (cell->live && cell->parent->solid <= 0 ) { // shake nothing inside a solid cell
			disp += get_rand(kick); // sum the displacement for each level
		} // but still apply the displacement from above
		cell->endN += disp;
		cell->xyz  += disp;
		cell->endC += disp;
	}
	// apply on the way down so high level displacment is passed down
	n = cell->kids;
	for (i=0; i<n; i++) shakeCell(cell->child[i],disp);
	//if (Data::shrink < NOISE) return;
	kick = 1.0 - Data::shrink;
	if (moltype == 3) { // 2D cells
		cell->xyz.z *= kick;
	} else { Vec was = cell->xyz;
		if (level==1) {
			cell->xyz  *= kick;
			cell->endN += cell->xyz-was;
			cell->endC += cell->xyz-was;
		}
	}
}

void spinsCell ( Cell *cell )
{
float	angle;
Vec	zero, axis;
int	i, m, n, id;
int	level = cell->level;
Data	*param = Data::model+cell->model;
int	moltype = param->moltype;
	n = cell->kids;
	if (cell->empty) return;
	angle = param->kicks[level];
	if (angle < 0.0) return; // -ve = no spin
	if (angle > NOISE) { // don't apply zero spin
		if (cell->live && cell->parent->solid <= 0 ) { // spin nothing inside a solid cell
			if (randf()<0.5) angle = -angle;
			cell->spin(angle);
		}
	}
	// high level rotation is passed down automatically in Cell::spin()
	DO(i,n) spinsCell(cell->child[i]);
}

void tinker ( float weight ) {
Cell	*world = Cell::world;
	DO(i,world->kids) // don't patch the world, just its offspring
	{ Cell	*child = world->child[i];
		if (child->live < 1) continue;
		if (child->parent->solid > 0) continue;
		if (child->empty) continue;
		patchCell(child,weight);
	}
}

void patchCell ( Cell *cell, float weight )
{
Data	*param = Data::model+cell->model;
int	level = cell->level;
int	dists = param->local[cell->level+1],
	bonds = param->chain[cell->level+1];
int	i,j,k, m,n, in;
int	fixes = 50;
float	dev, gev, w;
	n = cell->kids;
	in = 0;
	for (i=0; i<fixes && dists && bonds && n>3; i++) {
		in++;
		dev = gev = 0.0;
		// fix local distances
		m = (int)(5.0*randf());
		for (j=m; j<n; j+=5) {
			dev += fixdist(cell->child[j],weight*1.0);
		}
		// fix local angles
		m = (int)(5.0*randf());
		for (j=m; j<n; j+=5) {
			gev += fixgeom(cell->child[j],weight*1.0);
		}
		if (dev<0.1 && gev<0.1) break;
	}
	for (i=0; i<n; i++) patchCell(cell->child[i], weight);
}

float fixdist ( Cell *cell, float fix )
{
int	level = cell->level;
Data	*param = Data::model+cell->model;
float	bond = param->sizes[level]+param->bonds[level];
float	da1b1, dc0b2, dc0a2, da1b2, da2b2, db1a2;
Cell	*b2, *b1, *c0, *a1, *a2;
Cell	*parent = cell->parent, *jun = parent->starts, *sen = parent->finish;
float	f, d, dd = 0.0;
int	i, end=0;
	c0 = cell;
	b1 = c0->sis; a1 = c0->bro;
	if ( b1->level != level ) return 0.0;
	if ( a1->level != level ) return 0.0;
	if (b1->parent != parent) return 0.0;
	if (a1->parent != parent) return 0.0;
	b2 = b1->sis; a2 = a1->bro;
	if ( b2->level != level ) return 0.0;
	if ( a2->level != level ) return 0.0;
	if (b2->parent != parent) return 0.0;
	if (a2->parent != parent) return 0.0;
	// don't refine over free bonds (-ve prox) or zero length
	if (b2->prox.x<NOISE || b2->prox.y<NOISE || b2->prox.z<NOISE) return 0.0;
	if (b1->prox.x<NOISE || b1->prox.y<NOISE || b1->prox.z<NOISE) return 0.0;
	if (c0->prox.x<NOISE || c0->prox.y<NOISE || c0->prox.z<NOISE) return 0.0;
	if (a1->prox.x<NOISE || a1->prox.y<NOISE || a1->prox.z<NOISE) return 0.0;
	if (a2->prox.x<NOISE || a2->prox.y<NOISE || a2->prox.z<NOISE) return 0.0;
	// don't refine over cyclic ends
	if (b2==jun || b2==sen) end++;
	if (b1==jun || b1==sen) end++;
	if (c0==jun || c0==sen) end++;
	if (a1==jun || a1==sen) end++;
	if (a2==jun || a2==sen) end++;
	if (end > 1) return 0.0;
	// skip branched return connections
	if (c0==b2 || c0==a2 || b1==a1) return 0.0;
	// fix long distances (but not across parent level)
	da1b2 = c0->dist.x; // i-2-->i+1
	da2b2 = c0->dist.y; // i-2-->i+2
	db1a2 = c0->dist.z; // i-1-->i+2
	if (da2b2 > 0.0 && da2b2 < 999.0) {
		part2cells(a2,b2,da2b2,fix);
		d = vdif(a2->xyz,b2->xyz) - da2b2; dd += d*d;
	}
	if (da1b2 > 0.0 && da1b2 < 999.0) {
		part2cells(a1,b2,da1b2,fix);
		d = vdif(a1->xyz,b2->xyz) - da1b2; dd += d*d;
	}
	if (db1a2 > 0.0 && db1a2 < 999.0) {
		part2cells(b1,a2,db1a2,fix);
		d = vdif(b1->xyz,a2->xyz) - db1a2; dd += d*d;
	}
	// fix i--i+2 distances
	if (b1->geom.y < 999.0) {
		dc0b2 = b1->prox.y;
		part2cells(c0,b2,dc0b2,fix);
		d = vdif(c0->xyz,b2->xyz) - dc0b2; dd += d*d;
	}
	if (c0->geom.y < 999.0) {
		da1b1 = c0->prox.y;
		part2cells(a1,b1,da1b1,fix);
		d = vdif(a1->xyz,b1->xyz) - da1b1; dd += d*d;
	}
	if (a1->geom.y < 999.0) {
		dc0a2 = a1->prox.y;
		part2cells(c0,a2,dc0a2,fix);
		d = vdif(c0->xyz,a2->xyz) - dc0a2; dd += d*d;
	}
	// fix i--i+1 bonds
	part2cells(b1,b2,b1->prox.x,fix);
	part2cells(a2,a1,a1->prox.z,fix);
	if (randf()<0.5) {
		part2cells(a1,c0,c0->prox.z,fix);
	} else {
		part2cells(b1,c0,c0->prox.x,fix);
	}
	d = vdif(b1->xyz,b2->xyz) - b1->prox.x; dd += d*d;
	d = vdif(a2->xyz,a1->xyz) - a1->prox.z; dd += d*d;
	d = vdif(c0->xyz,a1->xyz) - c0->prox.z; dd += d*d;
	d = vdif(c0->xyz,b1->xyz) - c0->prox.x; dd += d*d;
	return dd;
}

float fixgeom ( Cell *cell, float fix )
{
int	level = cell->level;
Data	*param = Data::model+cell->model;
float	bond = param->sizes[level]+param->bonds[level];
float	tau1 = cell->geom.x, theta = cell->geom.y, tau2 = cell->geom.z;
float	dt1old,dthold,dt2old, dt1new,dthnew,dt2new;
Vec	x,y,z, mid, shift, xyz[5], abc[5];
Cell	*b2, *b1, *c0, *a1, *a2;
Cell	*parent = cell->parent, *jun = parent->starts, *sen = parent->finish;
float	dt, th, t1, t2, tt, good = 1.0;
float	d, dhi, dlo, tau;
Mat	mat, wat;
int	i, end=0;
	if (!cell->bond) return 0.0;
	if (tau1+theta+tau2 >10.0) return 0.0; // includes an unset angle (999)
	if (fix < NOISE) return 0.0;
	// set and check local cells: b2-b1-c0-a1-a2
	c0 = cell;
	b1 = c0->sis; a1 = c0->bro;
	if ( b1->level != level ) return 0.0;
	if ( a1->level != level ) return 0.0;
	if (b1->parent != parent) return 0.0;
	if (a1->parent != parent) return 0.0;
	b2 = b1->sis; a2 = a1->bro;
	if ( b2->level != level ) return 0.0;
	if ( a2->level != level ) return 0.0;
	if (b2->parent != parent) return 0.0;
	if (a2->parent != parent) return 0.0;
	// don't refine over free bonds (-ve prox) or zero length
	if (b2->prox.x<NOISE || b2->prox.y<NOISE || b2->prox.z<NOISE) return 0.0;
	if (b1->prox.x<NOISE || b1->prox.y<NOISE || b1->prox.z<NOISE) return 0.0;
	if (c0->prox.x<NOISE || c0->prox.y<NOISE || c0->prox.z<NOISE) return 0.0;
	if (a1->prox.x<NOISE || a1->prox.y<NOISE || a1->prox.z<NOISE) return 0.0;
	if (a2->prox.x<NOISE || a2->prox.y<NOISE || a2->prox.z<NOISE) return 0.0;
	// don't refine over cyclic ends
	if (b2==jun || b2==sen) end++;
	if (b1==jun || b1==sen) end++;
	if (c0==jun || c0==sen) end++;
	if (a1==jun || a1==sen) end++;
	if (a2==jun || a2==sen) end++;
	if (end > 1) return 0.0;
	d = vdif(b1->xyz,a1->xyz);
	dlo = c0->prox.y*(1.0-good);
	dhi = c0->prox.y*(1.0+good);
	if (d<dlo || d>dhi) return 99.9; // proximal distance too poor
	d = vdif(b2->xyz,a2->xyz);
	dlo = c0->dist.y*(1.0-good);
	dhi = c0->dist.y*(1.0+good);
	if (d<dlo || d>dhi) return 99.9; // distal distance too poor
	d = angle(b2->xyz,c0->xyz,a2->xyz)*180/PI;
	if (d < 30.0 || d > 150.0) return 0.0; // base triangle too thin
	t1 = torsion(b2->xyz,b1->xyz,c0->xyz,a1->xyz);
	if (t1>999.0) { Pt(Bad t1 in) Pi(cell->uid) NL
		Pv(b2->xyz) NL Pv(b1->xyz) NL Pv(c0->xyz) NL Pv(a1->xyz) NL exit(1);
	}
	th = angle(b1->xyz, c0->xyz, a1->xyz);
	t2 = torsion(a2->xyz,a1->xyz,c0->xyz,b1->xyz);
	if (t2>999.0) { Pt(Bad t2 in) Pi(cell->uid) NL
		Pv(a2->xyz) NL Pv(a1->xyz) NL Pv(c0->xyz) NL Pv(b1->xyz) NL exit(1);
	}
	dt1old = angdif(t1,tau1);
	dthold = angdif(th,theta);
	dt2old = angdif(t2,tau2);
	if (dt1old+dthold+dt2old > 3.0) return 99.9; // too far to fix
	xyz[0] = b2->xyz; xyz[1] = b1->xyz;
	xyz[2] = c0->xyz;
	xyz[3] = a1->xyz; xyz[4] = a2->xyz;
	// dial-up correct torsion angles
	t1 = torsion(xyz[3],xyz[2],xyz[1],xyz[0]);
	xyz[0].get_rot(xyz[2],xyz[1],tau1-t1);
	t2 = torsion(xyz[1],xyz[2],xyz[3],xyz[4]);
	xyz[4].get_rot(xyz[2],xyz[3],tau2-t2);
	t1 = torsion(xyz[0], xyz[1], xyz[2], xyz[3]);
	th = angle(xyz[1], xyz[2], xyz[3]);
	t2 = torsion(xyz[1], xyz[2], xyz[3], xyz[4]);
	dt1new = angdif(t1,tau1);
	dthnew = angdif(th,theta);
	dt2new = angdif(t2,tau2);
	// fit positions xyz[0,2,4] to b2,c0,a2
	x = xyz[4]- xyz[0];
	mid = xyz[4] & xyz[0];
	y = xyz[2] - mid;
	z = x^y;
	mat = Mat(x.norm(), y.norm(), z.norm()); // new unit basis vectors in mat
	wat = mat.get_inv();
	mid = (xyz[0] + xyz[2] + xyz[4])/3.0;	// new CoG
	for (i=0; i<5; i++) {  // get xyz from centre in terms of basis set coeficients (abc)
		xyz[i] -= mid;
		abc[i] = wat * xyz[i];
	}
	x = a2->xyz - b2->xyz;
	mid = a2->xyz & b2->xyz;
	y = c0->xyz - mid;
	z = x^y;
	x.setVec(); y.setVec(); z.setVec();		// old basis vectors
	mid = (b2->xyz + c0->xyz + a2->xyz)/3.0;	// old CoG
	for (i=0; i<5; i++) {	// reconstruct new positions with old basis vectors
		xyz[i] = mid;	// start at centre and sum basis vector components
		xyz[i] += x * abc[i].x;
		xyz[i] += y * abc[i].y;
		xyz[i] += z * abc[i].z;
	}

	t1 = torsion(xyz[0], xyz[1], xyz[2], xyz[3]);
	th = angle(xyz[1], xyz[2], xyz[3]);
	t2 = torsion(xyz[1], xyz[2], xyz[3], xyz[4]);
	dt1new = angdif(t1,tau1);
	dthnew = angdif(th,theta);
	dt2new = angdif(t2,tau2);
	// shift towards new positions (by fix) if all angles better
	if (dt1new>dt1old || dthnew>dthold || dt2new>dt2old) return 0.0;
	shift = (xyz[0]-b2->xyz)*fix; moveCell(b2,shift);
	shift = (xyz[1]-b1->xyz)*fix; moveCell(b1,shift);
	shift = (xyz[2]-c0->xyz)*fix; moveCell(c0,shift);
	shift = (xyz[3]-a1->xyz)*fix; moveCell(a1,shift);
	shift = (xyz[4]-a2->xyz)*fix; moveCell(a2,shift);
	t1 = torsion(b2->xyz,b1->xyz,c0->xyz,a1->xyz);
	th = angle(b1->xyz, c0->xyz, a1->xyz);
	t2 = torsion(a2->xyz,a1->xyz,c0->xyz,b1->xyz);
	t1 = angdif(t1,tau1);
	th = angdif(th,theta);
	t2 = angdif(t2,tau2);
	tt = t1 + th + t2;
	return tt;
}
/*
void patchPair ( Cells*, float );
void fixpair( Cells*, Cells*, Cells*, Cells* );
void fixup( Vec*, Vec*, Vec, int, float );

void patchPair ( Cells *cell, float weight )
{
int	level = cell->level;
int	i,j,k, m,n, in;
	n = cell->kids;
	if (n==0) return;
	for (i=0; i<n; i++) {
		patchPair(cell->child[i], weight);
	}
	if (cell->sort != 1) return;
	for (i=0; i<n-1; i++)
       	{ Cells *base2,*base1 = cell->child[i],
	       	*pair2,*pair1 = base1->link[0];
		if (!pair1) continue;
		if (pair1->parent != cell) continue;
		base2 = base1->bro;
		pair2 = base2->link[0];
		if (!pair2) continue;
		if (pair2->parent != cell) continue;
		fixpair(base1,base2,pair1,pair2);
	}
}

void fixpair ( Cells *b1, Cells *b2, Cells *p1, Cells *p2 )
{
	float	w = 0.1;
float	xy, xz, yz, dab, dcd;
Vec	ab, cd, a,b,c,d, x,y,z;
	vcopy(b1->xyz, &a); vcopy(b2->xyz, &b);
	vcopy(p1->xyz, &d); vcopy(p2->xyz, &c);
	dab = vdif(a,b); dcd = vdif(c,d);
	vave(a,b,&ab); vave(c,d,&cd);
	vsub(ab,cd,&z); vnorm(&z); // midpoint connection
	vsub(a,b,&x); vnorm(&x); // chain1 direction
	vsub(c,d,&y); vnorm(&y); // chain2 direction
	xy = vdot(x,y); xy *= xy;
	xz = vdot(x,z); xz *= xz;
	yz = vdot(z,y); yz *= yz;
	if (1){ //xy>xz && xy>yz) {
		fixup(&x,&y,z,1,w);
		vsub(a,x,&a); vadd(b,x,&b);
		vsub(c,y,&c); vadd(d,y,&d);
		vsub(a,b,&x); vnorm(&x); vmul(&x,dab*0.5);
		vsub(c,d,&y); vnorm(&y); vmul(&y,dcd*0.5);
		vadd(ab,x,&(b1->xyz)); vsub(ab,x,&(b2->xyz));		
		vsub(cd,y,&(p1->xyz)); vadd(cd,y,&(p2->xyz));		
	}
	//if (xz>xy && xz>yz) fixup(&x,&y,z,2);
	//if (yz>xy && yz>xz) fixup(&x,&y,z,3);
}

void fixup( Vec *x, Vec *y, Vec z, int fix, float w ) {
Vec	p, q;
	if (fix==1 || fix==3) { // xy or yz is worst
		vprod(z,*x,&p); // p is better y direction
		vsub(p,*y,y); vmul(y,w); // y = shift
	}
	if (fix==1 || fix==2) { // xy or xz is worst
		vprod(z,*y,&q); // q is better x direction
		vsub(q,*x,x); vmul(x,w); // x = shift
	}
}
*/
