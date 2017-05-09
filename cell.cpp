#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

Cell*   Cell::world;
Cell*   Cell::scene;
int     Cell::total;
int     Cell::allatoms;
Cell**  Cell::uid2cell;

float inEgg ( Vec, Vec, Vec, float );
float inEgg ( Vec, Seg, float );

// short temp names for Data::globals
char  *names;
int   *shape, *links, *chain, *split, *local, *align;
float *sizes, *bumps, *kicks, *keeps, *bonds, *repel, *rejel;
// common temp variables
int   total, depth, model, moltype, subtype;

void Cell::extend ( const int n, const int type, const int sort )
{ // extend an existing family up to n members
Cell    **tmpchild = new Cell*[kids];        // new list of pointers to children and ...
Cell    ***tmprank = new Cell**[kids];       // pointers to (pointers of) ranked children
	DO(i,kids) {
		tmpchild[i] = child[i];
                tmprank[i] = new Cell*[3];
               	DO(j,3)  tmprank[i][j] = rank[i][j];
	}
        child = new Cell*[n];        // list of pointers to children
        rank = new Cell**[n];        // pointers to (pointers of) ranked children
        DO(i,n) { Cell *ci;
                child[i] = new(Cell);
                rank[i] = new Cell*[3];
		ci = child[i];
               	DO(j,3) {
                       	rank[i][j] = ci;
                       	ci->ranks[j] = i;
               	}
		if (i<kids) {	// copy pointers to old children (not deleted)
			child[i] = tmpchild[i];
		} else {	// make new child
               		ci->parent = this;
               		ci->level = this->level+1;
               		ci->id = i;
			ci->type = type;
			ci->sort = sort;
			Cell::total++;		// count total cells created
			ci->uid = Cell::total;
			Cell::uid2cell[Cell::total] = ci;
		}
        }
        kids = n;
}

void Cell::spawn ( const int n, const int type, const int sort )
{ // create n children for the current cell (with type and sort)
        kids = n;
        child = new Cell*[kids];        // list of pointers to children
        rank = new Cell**[kids];        // pointers to (pointers of) ranked children
        DO(i,kids) { Cell *ci;
                child[i] = new(Cell);
                rank[i] = new Cell*[3];
		ci = child[i];
                ci->parent = this;
                ci->level = this->level+1;
                ci->id = i;
                DO(j,3) {
                        rank[i][j] = ci;
                        ci->ranks[j] = i;
                }
		ci->type = type;
		ci->sort = sort;
		Cell::total++;		// count total cells created
		ci->uid = Cell::total;
		Cell::uid2cell[Cell::total] = ci;
        }
}
void Cell::spawn ( const int n, const int type ) {
	spawn(n,type,0);
}
void Cell::spawn ( const int n ) {
	spawn(n,1,0); // default = spheres
}
void Cell::spawn () {
	spawn(1,1,0); // one sphere
}

void Cell::print () {
        for (int i=0; i<level; i++) PT
        Pi(uid) Pi(atom) Pi(btom) Pi(ctom) Pi(level) Pi(kids)
	Pi(type) Pi(sort) Pv(xyz) Pi(nbonds) Pi(nlinks) Pi(model) NL
	// cout << " ranks ="; for (int i=0; i<3; i++) cout << " " << ranks[i]; NL
        for (int i=0; i<kids; i++) child[i]->print();
}

void Cell::move () {
Vec	v = get_rand(1.0);
        xyz += v; endN += v; endC += v;
        for (int i=0; i<kids; i++) child[i]->move(v);
}

void Cell::move ( const float d ) {
Vec	v = get_rand(d);
        xyz += v; endN += v; endC += v;
        for (int i=0; i<kids; i++) child[i]->move(v);
}

void Cell::move ( const Vec v ) {
        xyz += v; endN += v; endC += v;
        for (int i=0; i<kids; i++) child[i]->move(v);
}

void Cell::rankXYZ ( const int n ) {
        for (int i=0; i<n; i++) { // make n passes
                for (int j=0; j<3; j++) sortXYZcell(j);
        }
}

///////// sorter

int Cell::sortXYZrank ( const int use )
{ // one pass bubble sort on the children of the current cell
float	ga, gb;
int	swaps = 0;
	for (int i=1; i<kids; i++)
	{ Cell *a = rank[i-1][use], *b = rank[i][use];
		if (a==0 || b==0) return 0;
		if (use==0) { ga = a->xyz.x; gb = b->xyz.x; }
		if (use==1) { ga = a->xyz.y; gb = b->xyz.y; }
		if (use==2) { ga = a->xyz.z; gb = b->xyz.z; }
		if (ga > gb) {
			rank[i][use] = a;
			rank[i-1][use]=b;
			a->ranks[use] = i;
			b->ranks[use] = i-1;
			swaps++;
		}
	}
	return swaps;
}

void Cell::sortXYZcell ( const int use )
{ // call sortXYZrank on the current cell and its offspring
	if (this==0) return;
	if (kids==0) return;
	if (empty) return;
	this->sortXYZrank(use);
	DO(i,kids) child[i]->sortXYZcell(use);
}

///////// utility

void moveCell ( Cell *cell, Vec disp, int shift )
{ // same as move() but keeps old style parameter list
int	ends = 1;
int	i, n = cell->kids;
	if (shift==0) return;
	if (disp.iszero()) return;
	if ((int)cell->endN.z==1234 || (int)cell->endC.z==1234) ends = 0;
	if (shift>0) {
		if (ends) cell->endN += disp;
		cell->xyz += disp;
		if (ends) cell->endC += disp;
	}
	if (shift<0) {
		if (ends) cell->endN -= disp;
		cell->xyz -= disp;
		if (ends) cell->endC -= disp;
	}
	if (n==0) return;
	for (i=0; i<n; i++) {
		moveCell(cell->child[i], disp, shift);
	}
}
void moveCell ( Cell *cell, Vec disp ) {
	moveCell(cell,disp,1);	// use positive shift
}

void Cell::spin ( const Vec zero, const Vec axis, const float theta )
{       // spin current cell around zero--axis vector (theta = radians)
Seg	line = Seg(zero,axis);
        xyz.set_rot(line,theta);
        if ((int)endN.z!=1234 && (int)endC.z!=1234) {
                endN.set_rot(line,theta);
                endC.set_rot(line,theta);
        }
        if (kids==0) return;
        DO(i,kids) child[i]->spin(zero,axis,theta);
}
void Cell::spin ( const Vec v, const float theta )
{       // spin current cell around centre--axis vector (theta = radians)
Vec	axis = v + xyz;
	spin(xyz,axis,theta);
}
void Cell::spin ( const Seg line, const float theta )
{       // spin current cell around the given line (theta = radians)
	spin(line.A,line.B,theta);
}
void Cell::spin ( const float theta )
{       // spin current cell around a random axis (theta = radians)
Vec	axis = xyz + get_rand();
	spin(xyz,axis,theta);
}

void spinCell ( Cell *cell, const Vec zero, const Vec axis, const float theta )
{       // spin current cell around zero--axis vector (old style param list)
	cell->spin(zero,axis,theta);
}

/*
void turnCell ( Cells *cell, Vec zero, Mat *rot )
{
Vec     new;
int     i, n = cell->kids;
        vsub(cell->xyz,zero, &new);
        MmulV(rot,new, &new);
        vadd(new,zero, &(cell->xyz));
        if ((int)cell->endN.z!=1234 && (int)cell->endC.z!=1234) {
                vsub(cell->endN,zero,&new);
                MmulV(rot,new, &new);
                vadd(new,zero,&(cell->endN));
                vsub(cell->endC,zero,&new);
                MmulV(rot,new, &new);
                vadd(new,zero,&(cell->endC));
        }
        if (n==0) return;
        for (i=0; i<n; i++) turnCell(cell->child[i],zero,rot);
}
*/

void part2cells ( Cell *a, Cell *b, float dist, float kick )
{
Seg	ends;
	if (dist < 0.0) return;
	if (dist > 9999.0) return;
	if (a==0 || b==0) { Pt(Bad data) Pi(a->uid) Pi(b->uid) NL return; }
	if (a->level != b->level) return;
	ends.A = a->xyz; ends.B = b->xyz;
	ends.separate(dist,kick);
	a->move(ends.A - a->xyz);
	b->move(ends.B - b->xyz);
}	
void part2cells ( Cell *a, Cell *b, float dist ) {
	part2cells(a,b,dist,1.0);	// move to full separation
}

Vec Cell::getWcent ()
// get the centroid of the children weighted by their number of kids
{
Vec	cog, add;
float	wkids = 0.0;
	if (kids < 1) return xyz;
	if (kids < 2) return child[0]->xyz;
	cog.zero();
	DO(i,kids)
        { float w = 1.0+(float)(child[i]->kids); // weight into sum by kids
		cog += child[i]->xyz * w;
               	wkids += w;
        }
        cog /= wkids;            // find weighted centre of children
	return cog;
}

void Cell::setWcent () {
// not substituted directly to avoid cell->xyz=0
	this->xyz = this->getWcent();
}

///////// bumper

float touch ( Cell *a, Cell *b )
{ // closest approach between two cell bump surfaces (of the same shape type)
  // +ve = separation, -ve = penetration depth
float	ra, rb, da, db, d, rab;
int	ta, tb, tab;
Vec	pa, pb;
Seg	s;
	if (a->type > b->type) { Cell *c = a; a = b; b = c; } // swap
	ta = a->type; tb = b->type;
	da = Data::model[a->model].bumps[a->level]; ra = da*0.5;
	db = Data::model[b->model].bumps[b->level]; rb = db*0.5;
	pa = a->xyz; pb = b->xyz;
	rab = ra + rb;
	tab = ta * tb;
	if (tab>1) s = Seg(b->endN,b->endC);
	switch (tab) { // smallest object type first
		case 0:	// virtual spheres
			return (pa|pb)-rab; // 9999.9;
		case 1:	// spheres = centre distance
			return (pa|pb)-rab;
		case 2:	// sphere+tube (closest approach to line segment or ends)
			if (pa.vec_in_seg(s)) return pa.vec_to_line(s)-rab;
			return fmin(pa|s.A,pa|s.B) - rab;
		case 3:	// sphere+ellipsoid (in keeper.cpp)
			d = inEgg(pa,s,db) - ra;
			if (d < 0.0) return d;
			if (d > rab) return d;
			return vec_to_egg(pa,s,db) - ra;
		case 4:	// tubes = closest approach of 2 line segments
			return seg_to_seg(Seg(a->endN,a->endC),Seg(b->endN,b->endC))-rab;
		case 6:	// tube+ellipsoid (in bumper.cpp)
			return tube_to_egg(a,b);
		case 9: // ellipsoid (in bumper.cpp)
			return egg_to_egg(a,b);
	}
}

int bumpex ( Cell *a, Cell *b ) {
// the children of the cell <a> and <b> are checked for inter-family bumps
int	i, j, level = a->level, kidlev, n = 0;
float	d, bumpa,bumpb, bump, over = 1.101;
float	kicka, kickb, hard, soft, kick;
float	strength = 0.1;
Vec	axis;
	if (a->solid > 0 || b->solid > 0 ) return 0;
	if (a->empty > 0 || b->empty > 0 ) return 0;
	if (a->kids == 0 || b->kids == 0 ) return 0;	// no children to bump
	if (a->level != b->level) return 0; // different level
	kidlev = level+1;
	kicka = Data::model[a->model].repel[kidlev],
	kickb = Data::model[b->model].repel[kidlev],
	hard = 0.5*(kicka+kickb);
	kicka = Data::model[a->model].rejel[kidlev],
	kickb = Data::model[b->model].rejel[kidlev],
	soft = 0.5*(kicka+kickb);
	if (kidlev==Data::depth) kick = hard; else kick = soft;
	axis = b->xyz - a->xyz;
	axis.setVec(soft);	// soft length vector from a to b (NB has to be set at atom level)
	DO(i,a->kids) { Cell* ai = a->child[i];
		if (ai->empty) continue;
		DO(j,b->kids) { Cell* bj = b->child[j];
			if (bj->empty) continue;
			if (exempt(ai,bj)) continue;
			bump = touch(ai,bj);
			if (bump > -NOISE) continue;
			d = (ai->xyz|bj->xyz)-bump*over; // d = target gap (NB clash has -ve bump)
			moveCell(ai,axis,-1);	// nudge ai towards a
			moveCell(bj,axis, 1);	// nudge bj towards b
			part2cells(ai,bj,d,-kick*strength); // -kick = repel only
                        ai->bump = LIVE; bj->bump = LIVE;
                        ai->hit = bj; bj->hit = ai;
			n++;
		}
	}
	return n;
}

int Cell::bumpin () {
// the children of the current cell <this> are checked for intra-family bumps
Bumps	*list;
int	i, j, kidlev, in, m, n = 0;
float	d, bump, boot, kick, over = 1.01;
float	hard, soft;
int	weight = 0;
	if (this->bump > 0) { // count-down refractory period and if over, clear hit
	// the plan is to use the period of bumping to activate MD/refinement in potter() 
		this->bump--;
		if (this->bump==0) this->hit = 0;
	}
	if (solid > 0 ) return 0;
	if (empty > 0 ) return 0;
	if (kids == 0 ) return 0;	// no children to bump
	list = new Bumps[HOLD*kids];
	kidlev = level+1;
	hard = Data::model[model].repel[kidlev]; // hard (for hard-shell bump)
	soft = Data::model[model].rejel[kidlev]; // soft (for jelly bumping)
	depth = Data::depth;
	if (kidlev==depth) {
		kick = hard; // atomic level has hard-shell bump
	} else {
		kick = soft; // higher levels have jelly bumping modified by kids if...
		if (kick > NOISE) weight = 1; // +ve = weight kick by kids (as Gauss(m) soft<-->hard)
	}
	if (kick < 0.0) kick = -kick;
	in = getBumpin(this,list); // list of bumping pairs (by NxN or sort for big N) 
	n = 0;
	for (i=0; i<in; i++) { Cell *a = list[i].a, *b = list[i].b;
		bump = touch(a,b);
		if (bump > -NOISE) continue;
		if (kidlev<depth) m = bumpex(a,b); else m = 0;	// bumping a+b children parted in bumpex() 
		if (exempt(a,b)) continue;	// exempt parents (exempt atoms skipped in getBumpin()) 
		// the pair (a,b) are bumping so repel more with more bumping children (m)
		//	unless weight=0 then just use unmodified <soft> value
		if (weight) { // Gaussian switch from soft to hard with increasing <m>
			d = (float)m; d = exp(-d*d);
			boot = d*soft + (1.0-d)*hard;
			boot *= kick;
		} else { boot = kick; }
		d = (a->xyz|b->xyz)-bump*over; // clash = -ve bump
		part2cells(a,b,d,-boot); // -kick = repel only
		if (kidlev<depth) {
               		a->bump = LIVE; b->bump = LIVE;
                        // a->hit = b; b->hit = a; // replace current interaction? trapped in exempt()?
               		if (!a->hit) a->hit = b;
			if (!b->hit) b->hit = a;
		}
		n++;
	}
	delete [] list;
}

void Cell::bumps () {
	this->bumpin();
	for (int i=0; i<kids; i++) child[i]->bumps();
}

///////// keeper

void Cell::shifter ( const float couple )
{ // shift child centroid to cell centre by couple and cell to centroid by 1-couple
Vec	shift, s;
float	wkids = 0.0;
	if (kids==0) return;
	shift = xyz - this->getWcent();		// shift for children towards cell centre
	s = shift * couple;			// scale by coupling constant
	DO(i,kids) moveCell(child[i],s,1);	// move children (and all their offspring) home
	s = shift * (1.0 - couple);		// scale by 1-couple
	xyz -= s; endN -= s; endC -= s;		// shift shell (parent and poles) towards children
}
void Cell::group () {
	this->shifter( 0.5 );
	for (int i=0; i<kids; i++) child[i]->group();
}

bool kidsLink ( Cell *A, Cell *B ) {
// returns TRUE if any of the children of A and B are linked
	if (A->kids==0) return FALSE;
	if (B->kids==0) return FALSE;
	// check if chlid level has any links
	if (A->child[0]->link==0) return FALSE;
	if (B->child[0]->link==0) return FALSE;
	// check if A kids link to B kids
	DO(i,A->kids) { Cell *c = A->child[i];
		if (c->nlinks==0) continue;
		DO(j,c->nlinks) {
			if (c->link[j].to==0) continue;
			if (c->link[j].to->parent == B) return TRUE;
		}
	}
	// check if B kids link to A kids
	DO(i,B->kids) { Cell *c = B->child[i];
		if (c->nlinks==0) continue;
		DO(j,c->nlinks) {
			if (c->link[j].to==0) continue;
			if (c->link[j].to->parent == A) return TRUE;
		}
	}
	return FALSE;
}

bool kidsLink ( Cell *A, int ia, Cell *B ) {
// returns TRUE if child ia of A and any child of B are linked
	if (A->kids==0) return FALSE;
	if (B->kids==0) return FALSE;
	// check if chlid level has any links
	if (A->child[0]->link==0) return FALSE;
	if (B->child[0]->link==0) return FALSE;
	// check if A kids link to B kids
	DO(i,1) { Cell *c = A->child[ia];
		if (c->nlinks==0) continue;
		DO(j,c->nlinks) {
			if (c->link[j].to==0) continue;
			if (c->link[j].to->parent == B) return TRUE;
		}
	}
	// check if any B kids link to Aa kid
	DO(i,B->kids) { Cell *c = B->child[i];
		if (c->nlinks==0) continue;
		DO(j,c->nlinks) {
			if (c->link[j].to==0) continue;
			if (c->link[j].to == A->child[ia]) return TRUE;
		}
	}
	return FALSE;
}


bool kidsLink ( Cell *A, int ia, Cell *B, int ib ) {
// returns TRUE if child ia of A and child ib of B are linked
	if (A->kids==0) return FALSE;
	if (B->kids==0) return FALSE;
	// check if chlid level has any links
	if (A->child[0]->link==0) return FALSE;
	if (B->child[0]->link==0) return FALSE;
	// check if Aa kid links to any B kids
	DO(i,1) { Cell *c = A->child[ia];
		if (c->nlinks==0) continue;
		DO(j,c->nlinks) {
			if (c->link[j].to==0) continue;
			if (c->link[j].to == B->child[ib]) return TRUE;
		}
	}
	// check if Bb kid links to Aa kid
	DO(i,1) { Cell *c = B->child[ib];
		if (c->nlinks==0) continue;
		DO(j,c->nlinks) {
			if (c->link[j].to==0) continue;
			if (c->link[j].to == A->child[ia]) return TRUE;
		}
	}
	return FALSE;
}

void putall ( float s )
{ // dump everything (grouped by level with scale = s)
FILE	*out = fopen("temp.pdb","w");
	DO1(j,Data::depth) {
		DO1(i,Cell::total) { Cell *c = Cell::uid2cell[i]; int lev = c->level;
			if (lev != j) continue;
			fprintf(out,"ATOM%7d  CA  GLY %c%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", c->atom, 'A'+lev, c->btom,
                        	s*c->xyz.x, s*c->xyz.y, s*c->xyz.z, (float)c->type, (float)c->sort);
		}
		fprintf(out,"TER ATOM\n");
	}
}

void putall ()
{ // dump everything (grouped by level with s = Data::scaleout)
	putall(Data::scaleout);
}

void putpdb ( Cell *top )
{ // dump each chain (atomic level only) into "temp.pdb"
	putpdb("temp.pdb",top,Data::scaleout);
}

void putpdb ( Cell *top, float s )
{ // dump each chain (atomic level only) into "temp.pdb"
	putpdb("temp.pdb",top,s);
}

void putpdb ( char *file, Cell *top, float s )
{ // dump each chain (atomic level only) into <file>
FILE	*out = fopen(file,"w");
	putpdb(out,top,s);
	fclose(out);
}

void putpdb ( char *file, float s )
{ // dump each chain (atomic level only) into <file>
FILE	*out = fopen(file,"w");
	putpdb(out,Cell::world,s);
	fclose(out);
}

void putpdb ( char *file )
{ // dump each chain (atomic level only) into <file>
FILE	*out = fopen(file,"w");
	putpdb(out,Cell::world,Data::scaleout);
	fclose(out);
}

void putpdb ( FILE *out, Cell *top, float s )
{ // dump each chain (atomic level only) to <out>
float	sec;
int	m = 0, n = 0;
	depth = Data::depth;
	DO1(i,Cell::total)
	{ Cell	*p, *c = Cell::uid2cell[i];
		if (c->level != depth) continue;
		p = c;
		for (int j=depth; j>top->level; j--) p = p->parent; 
		if (p != top) continue;
		sec = 0.0;
		if (c->parent->sort==1) sec = 2.0;
		if (c->parent->sort==2) sec = 1.0;
		if (m && c->btom==1) {
			n++;
			fprintf(out,"TER\n");
		}
		fprintf(out,"ATOM%7d  CA  GLY %c%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", ++m, 'A'+n, c->resn,
                       	s*c->xyz.x, s*c->xyz.y, s*c->xyz.z, (float)c->type, sec);
	}
}
