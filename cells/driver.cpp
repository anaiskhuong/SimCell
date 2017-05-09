#include "main/util.hpp"
#include "main/geom.hpp"
#include "main/cell.hpp"
#include "main/data.hpp"

int run;
int history  = 0;
int silent   = 0;

// model control feature
int bumping = 1; // 0 = bumping off, otherwise on
int bolting = 1; // green cells into fast escape mode  (used as <bolts>)
int fleeing = 1; // green cells become more persistent (used as <flees>)
int hooking = 1; // 1=all, 2=RR+GG, 3=RR, 4=GG, 5=RG   (used as <hooks>)
int turning = 0; // reverse lead direction on bumping  (used as <turns> set by -ve bumping value)
int memory  = 0; // memory of recent encounters
// apply to cell type (model) <turns>-1 and if over <nmodels>, to all
// motion control parameters
float	Decay  = 0.95;		// speed decay after bumping
float	Slows  = 0.99;		// speed decay after hooking
float	Unhook = 0.001;		// 0.001; // chance to unhook
int	Flow = 0;		// switch on/off leader contact inhibition
float	Wild = 0.0;		// fast speed for mixing
float	Pull = 0.0;		// directed kick to leader (speed)
float	Fast = 0.0;		// factor (0..1) to modify <Speed> in panic mode
float	Duck = 0.0;		// chance of reassigning the leader
float	Dive = 0.0;		// factor (0..1) to modify <Dodge> in panic mode
float	Speed, Dodge;		// temp values of <Pull> and <Duck> 
int	Split  = 0;		// split fields in half
int	Finger = 0;		// extend a finger into each half
int	Freeze = 0;		// freeze the motion of cell type 0 (set by -Pull)
int Fuzzy = 0;

float xfinger=4.0, yfinger=20.0;
float yposition = 20.0;
int poke = 16;

// global values used by looker()
int	*domin;
int	print, cycles, delays, hold, length, start;
FILE	*pdb, *msd, *rdf, *dom, *coord;
int	*countMSD, maxMSDrun = 1000;
float	*sumdd;

int bolts, flees, hooks, turns;

void busyBody ( Cell* );
void busyBead ( Cell* );
void setLinks ( Cell* );
void memories ( Cell* );
int  newLeader ( Cell*, int );



void looker () { // called on same thread as helper
	Cell	*world = Cell::world;
	int	n = world->kids;
	float   edge = 0.5*Data::model[0].sizes[0]; // world's radius
	float	h[999];
	int	p[999], q[999];
	Vec	m[999], last;
	int	in = 1; // counts close pairs within range
	float	len = 0.0;
	char	file[22];
	FILE	*mid;
	int	frame = Data::frame;
	sprintf(file,"split%d.dat", frame/1000);
	mid = fopen(file,"a");
	sleep(5);
	if (Cell::world->done == 0) return; // wait
	if (Data::frame < 300) return;
	m[0].x = m[0].z = 0; m[0].y = h[0] = -edge;
	DO(i,n) { Cell *ci = world->child[i]; float d, dmin = 99999.9;
		DO(j,n) { Cell *cj = world->child[j];
			if (ci->model == cj->model) continue;
			d = ci->xyz|cj->xyz;
			if (d > dmin) continue;
			if (d > 5.0) continue;
			dmin = d;
			ci->vdata[0] = ci->xyz;
			ci->vdata[1] = cj->xyz;
		}
		if (dmin > 99999.0) continue;
		m[in] = ci->vdata[0] & ci->vdata[1];
		h[in] = m[in].y;
		in++;
	}
	m[in].x = m[in].z = 0; m[in].y = h[in] = edge;
	in++;
	//sort(h,p,in); // sort mid-points on y
	// follow down
	DO(i,in) h[i] = -1; // flag for not used (else 1)
	n = 0; p[0] = n; h[n] = 1; //  start at top
	DO(i,in) { int next = -1; float d, dmin = 999.9;
		if (i==0) continue;
		DO(j,in) { // find closest point to current (n)
			if (h[j] > 0) continue; // used
			d = m[n] | m[j];
			if (d > dmin) continue;
			dmin = d;
			next = j;
		}
		n = next; p[i] = n; h[n] = 1;
	}
	// follow up
	DO(i,in) h[i] = -1; // flag for not used (else 1)
	n = in-1; q[0] = n; h[n] = 1; //  start at bottom
	DO(i,in) { int next = -1; float d, dmin = 999.9;
		if (i==0) continue;
		DO(j,in) { // find closest point to current (n)
			if (h[j] > 0) continue; // used
			d = m[n] | m[j];
			if (d > dmin) continue;
			dmin = d;
			next = j;
		}
		n = next; q[i] = n; h[n] = 1;
	}
	last = m[0];
	for(int i=0;i<(in-1);i++)
	{ int	j = p[i], k = q[in-i-1];
		float fin = (float)(in-1), wq = (float)i/fin, wp = (float)(in-i-1)/fin;
		Vec	mpq = m[j]*wp + m[k]*wq;
		fprintf(mid,"%f %f\n", -mpq.x,mpq.y);
		if (i>0) { len += mpq|last; last = mpq; }
	}
	fprintf(mid,"\n");
	fclose(mid);
	Pi(frame) Pr(len) NL
}

void driver ()
{
    Cell	*world = Cell::world;
    Data	*model = Data::model+world->model;
    int	n = world->kids;
    int	split;
    char filename[256];
    run = Data::frame;
    if (run%1000==0) { Pi(run) NL }
    if (run < 1) return;
    if (world->idata[0]==0) { FILE *dat;
        world->idata[0] = 1; // flag started
        print = 0;
        domin = new int[n+1];
        DO(i,n) domin[i] = 0;
        countMSD = new int[maxMSDrun];
        DO(i,maxMSDrun) countMSD[i] = 0;
        sumdd = new float[maxMSDrun];
        DO(i,maxMSDrun) sumdd[i] = 0.0;
        world->len = model->sizes[0];
        dat = fopen("cell.param","r");
        if (dat) {
            printf("Parameter file read\n");
            fscanf(dat,"%d %d %d %d", &cycles, &delays, &hold, &length);
            //			delays += n/5;
            Pi(cycles) Pi(delays) Pi(hold)Pi(length) NL
                    fscanf(dat,"%d %d %d %d", &bumping, &hooking, &bolting, &fleeing);
            if (bumping < 0) { turning = -bumping; bumping = 1; }
            if (hooking < 0) { hooking = -hooking; memory = 1; }
            Pi(bumping) Pi(hooking) Pi(bolting) Pi(fleeing) Pi(turning) Pi(memory) NL
                    fscanf(dat,"%f %f %f %f %f %f %f", &Wild, &Pull, &Fast, &Duck, &Dive, &Unhook, &Decay);
            // <Wild> is a faster speed for initial randomisation
            // <Fast> is a factor to modify <Speed> (default <Pull>) big = faster
            // <Dive> is a factor to modify <Dodge> (default <Duck>) big = straighter
            // <Unhook> id the chance to break a link (not active)
            // >Decay> rate of memory decay
            if (Pull < 0.0) { Freeze = 1; Pull = -Pull; }
            if (Decay < 0.0) { Flow = 1; Decay = -Decay; }
            if (Wild < 0.0) {
                Split = 1; Wild = -Wild;
                if ((int)Wild == 10) { Finger = 1; Wild = 0.01; }
            }
            if (Finger) Split = 1;
            Pr(Wild) Pr(Pull) Pr(Fast) Pr(Duck) Pr(Dive) Pr(Unhook) Pi(Decay) NL
                    Pi(Flow) Pi(Split) Pi(Finger) Pr(Fuzzy) NL
                    fclose(dat);
            DO(i,n)	// loop over cells
            { Cell *cell = world->child[i];
                int	it;
                cell->fdata[0] = 1.0;		// fdata[0] = current <speed> factor
                cell->fdata[1] = 1.0;		// fdata[1] = current <dodge> factor
                cell->fdata[2] = 0.0;		// fdata[2] = max bead score
                cell->fdata[3] = 0.0;		// fdata[3] = min bead score
                cell->fdata[4] = 0.0;		// fdata[4] = max-min bead
                cell->busy = it = (int)(drand48()*(float)(cell->kids)); // busy=<leader>
                cell->child[it]->busy = 1;	// flag 1/0 for leader at bead level
                cell->idata[1] = 0;		// idata[1] = not used
                cell->idata[2] = 0;		// idata[2] = time in excited state
                DO(j,cell->kids) {	// loop over beads
                    // busy = leader
                    // idata[]: counters  0 = leader-link, 1..4 No.links to cell types
                    // fdata[]: memories  0 = leader-link, 1..4 No.links to cell types
                    DO(k,NDATA) {
                        cell->child[j]->idata[k] = 0;	// zero linkage counts
                        cell->child[j]->fdata[k] = 0.0;	// zero linkage memory
                    }
                    cell->child[j]->colour = -1;	// default col = cell
                    cell->child[j]->nlinks = 0;	// counts current links (not slots)
                }
            }
        } else {
            Pt(No cell.param file found) NL exit(1);
        }
        bolts = flees = hooks = turns = 0;	// switch off all behaviour modes
        Speed = Wild; // <Speed> is reset to default <Pull> below
        Dodge = Duck; // <Dodge> does not get changed
		Cell::world->done = 1; // go
    }
	
    if (run==delays/4) {
        //Split = 0;		// stop cell separation
        Speed = Pull;		// set speed back to normal rate
        if (Split) Speed = 0.0; // set speed to zero
        printf("wild time stopped at %d\n", run);
    }
    if (run==delays) {
        printf("settling time stopped at %d\n", run);
        Speed = Pull;		// set speed back to normal rate
        if (Split) Speed = 0.0;	// set speed to zero
        start = delays;
        print = 2;	// switch on start sample (2 = with pdb)
        if (Finger) print = 0;
        printf("start sample kept from %d to %d\n", delays, delays+hold);
        DO(i,n) { Cell *ci = world->child[i];
                  ci->vdata[0] = ci->xyz;	// record starting position
                }
        pdb = fopen("start.pdb","w");
        msd = fopen("start.msd","w");
        rdf = fopen("start.rdf","w");
        dom = fopen("start.dom","w");
    }
    if (Split) { float reflect = -3.0;
        DO(i,n) { Cell *ci = world->child[i];
                  if (Finger) {
                      if (i<poke) continue;
                      if (i>n/2 && i<n/2+poke) continue;
                      if (run < delays+hold+100)
                      {
                          float d = 1.0/(1.0+ci->xyz.x);
                          Vec	v = Vec(0,d*ci->xyz.y,0);
                          if (ci->model==1
                                  && ci->xyz.y<(yposition+yfinger/2-(yfinger/2)/xfinger*ci->xyz.x)
                                  && ci->xyz.y>(yposition-yfinger/2+(yfinger/2)/xfinger*ci->xyz.x))
                          {
                              ci->move(-v);
                          }
                          if (ci->model==0
                                  && ci->xyz.y<(-(yposition-yfinger/2)+(yfinger/2)/xfinger*ci->xyz.x)
                                  && ci->xyz.y>(-(yposition+yfinger/2)-(yfinger/2)/xfinger*ci->xyz.x))
                          {
                              ci->move( v);
                          }
                      }
                  }
                  if (run < delays+hold) {
                      if (ci->model) {
                          if (ci->xyz.x < 0.0) ci->xyz.x *= reflect;
                      } else {
                          if (ci->xyz.x > 0.0) ci->xyz.x *= reflect;
                      }
                  }
                }
    }
    if (Finger && run<delays+hold && run%5==0)
    { int	k = n/2+1;
        DO(i,poke)
        {
            int j = i+k; float x = xfinger, y = yfinger/2;
            // green
            float alpha = randf();
            world->child[i]->xyz.x = alpha*x;
            world->child[i]->xyz.y = yposition + (randf()*2-1)*(y-y*alpha);
            // red
//            world->child[j]->xyz.x = -x*(float)i+2.0;
//            world->child[j]->xyz.y = -y-randf()*yposition;
            world->child[j]->xyz.x = -alpha*x;
            world->child[j]->xyz.y = -(yposition + (randf()*2-1)*(y-y*alpha));
        }
        DO(i,poke) DO(j,8) {
            world->child[i]->child[j]->xyz.z = 0.0;
            world->child[i+k]->child[j]->xyz.z = 0.0;
        }
    }
    if (Fuzzy && run==delays+hold)
    {
        int swap, border, tmpmodel;
        int	p[999], q[999];
        float swapprob=0.3;
        DO(k,2){
            border = 0;
            DO(i,n) {
                Cell *ci = world->child[i];
                float d, dmin = 99999.9;
                {
                    DO(j,n) {
                        Cell *cj = world->child[j];
                        if (ci->model == cj->model) continue;
                        d = ci->xyz|cj->xyz;
                        if (d > dmin) continue;
                        if (d > 3.0) continue;
                        dmin = d;
                        p[border] = ci->id;
                        q[border] = cj->id;
                        border++;
                    }
                }
            }
            if (border)
            {
                swap=0;
                DO(i,border) {
                    if (randf()<(swapprob)){///(float)(k+1))){
                        Cell *ci = world->child[p[i]];
                        Cell *cj = world->child[q[i]];
                        tmpmodel = cj->model;
                        ci->model = cj->model;
                        DO (l,ci->kids){
                                    Cell *cl = ci->child[l];
                                    cl->model = ci->model;
                        }
                        swap++;
                    }
                }
            }
        }
        Fuzzy=0;
    }

    if (run==delays+hold) {
        Finger = 0;
        print = -1;	// switch off start sample (-1 = close files)
        printf("full behaviour started at %d\n", run);
        Speed = Pull;		// set speed back to normal rate
        hooks = hooking;
        flees = fleeing;
        bolts = bolting;
        turns = turning; 
    }
    if (run==delays+hold+cycles) {
        start = delays+hold+cycles;
        print = 2;	// switch on final sample (2 = with pdb)
        printf("final sample kept from %d to %d\n", run,run+hold);
        DO(i,n) { Cell *ci = world->child[i];
                  ci->vdata[0] = ci->xyz;	// record starting position
                }
        pdb = fopen("final.pdb","w");
        msd = fopen("final.msd","w");
        rdf = fopen("final.rdf","w");
        dom = fopen("final.dom","w");
    }
    if (run > delays+hold+cycles+hold) { FILE *done;
        print = -1;	// switch off final sample (-1 = close files)
        printf("Stopping the simulation\n");
        done = fopen("done","w");
        fprintf(done,"done\n");
        fclose(done);
        sleep(2); // wait for looker() to finish
        exit(1);
    }
    DO(i,n) { int id = world->rank[i][2]->id; // loop in order of rank in Z (random)
              busyBody(world->child[id]);
            }
    //////////////////

    if ( print == -1 && (run-delays-hold)%length == 0 ) {
        int linkDiffCell, linkedCells = 0;
        sprintf(filename,"cells_coord-%08d.txt",(run-delays-hold));
        coord = fopen(filename,"w");
        DO(i,n)
        {
            Cell *ci = world->child[i];
            linkDiffCell = 0;
            DO (j,ci->kids){
                Cell *cj = ci->child[j];
                if (cj->nlinks>0) {
                    linkedCells = 1;
                    DO (k,cj->nlinks){
                        if (ci->model != cj->link[k].to->model){linkDiffCell = 1;}
                    }
                }
            }
            fprintf(coord, "%2.2f %2.2f %2.2f %d\n", ci->xyz.x, ci->xyz.y, ci->xyz.z, ci->model);
            linkedCells = 0;
        }
        fclose(coord);
    }

    /////////////////
}

void busyBody ( Cell *cell ) {
    // animate the main cell body
    float	edge = 0.5*Data::model[0].sizes[0]; // world's radius
    float	speed  = cell->fdata[0],
            dodge  = cell->fdata[1];
    int	leader = cell->busy,
            n = cell->kids,
            id = cell->id,
            m = n-1;
    Cell	*lead = cell->child[leader];
    int	lea, der;
    Vec	kick;
    float	panic, escape, change;
    cell->xyz.z *= 0.9;
    // <Decay> = (0..1): 0.1 = fast, 0.9 = slow return to defalut
    if (speed < 1.0) speed /= Decay; 	// speed up towards normal
    if (speed > 1.0) speed *= Decay;	// slow speed towards normal
    // <speed> is a factor (0..1..)  that modifies Speed
    if (dodge < 1.0) dodge /= Decay; 	// revert to normal persistence
    if (dodge > 1.0) dodge *= Decay;	// revert to normal persistence
    // <speed> is a factor (0..1..)  that modifies Speed
    panic = 1.0;
    if (memory) { float f = cell->fdata[4];
        // use leader bead value
        if (lead->fdata[4] > 0.0) panic = 1.0+lead->fdata[4];
        // use cell average value
        //if (f > 0.0) panic = Fast*(1.0-exp(-f*f));
        //if (f < 0.0) panic = exp(-f*f);
    }
    //	if (Freeze==1 && cell->model==0) speed = 0.0; // freeze cell type-0 motion
    if (Freeze==1 && cell->model==0) {
        speed = 0.0; // freeze cell type-0 motion
        cell->xyz *= 0.99;
    }
    escape = Speed*speed;
    if (bolts) {
        escape *= panic;
    }
    change = Dodge/dodge;
    if (flees) {
        change /= panic;
    }
    lead->busy = 0;			// clear flag on old leader
    if (lead->xyz.len() > edge*0.9) {	// near the edge of the world
        lead->fdata[0] = 50.0;		// frighten the leader away from the edge
        leader = newLeader(cell,999);	// applies to all cell types (999=edge flag)
    } else {				// otherwise take a chance on a new leader
        if (randf()<change)
        {
            leader = newLeader(cell,0); // big dodge = fewer changes
        }
    }
    lead = cell->child[leader];
    lead->busy = 1;			// set flag for new leader
    cell->bump = 0;				// clear bump flag
    kick = lead->xyz - cell->xyz;		// kick = leader direction
    kick.setVec(escape);			// set kick size to modified <Speed>
    lead->xyz += kick;			// kick the leader
    cell->xyz += kick;			// and body
    kick *= 0.5;
    lea = leader-1; if (lea<0) lea=m;
    cell->child[lea]->xyz += kick;		// half kick to leader's left
    der = leader+1; if (der>m) der=0;
    cell->child[der]->xyz += kick;		// half kick to leader's right
    cell->fdata[0] = speed;
    cell->fdata[1] = dodge;
    cell->busy = leader;
    DO(i,n) setLinks(cell->child[i]);		// add new links and break long links
    if (memory) {
        DO(i,n) memories(cell->child[i]);	// update the linkage memory and colour
    }
    if (memory)
    { float dd, msd, temp[8], spread = 0.01;	// spread 1/10 to cyclic adjacent beads
        int	blue = 0, blac = 0;
        // smooth memories
        DO(j,5) { // for each memory type
            DO(i,n)	// share memories with neighbours
            { Cell	*bi, *ci, *di;
                float last, curr, next;
                if (i == 0) bi = cell->child[n-1]; else bi = cell->child[i-1];
                ci = cell->child[i];
                if (i< n-1) di = cell->child[i+1]; else di = cell->child[0];
                last = bi->fdata[j];
                curr = ci->fdata[j];
                next = di->fdata[j];
                temp[i] = curr*(1.0-spread)+spread*0.5*(last+next);
            }
            DO(i,n) cell->child[i]->fdata[j] = temp[i]; // copy in smoothed values
        }
        DO(i,n) { Cell *bead = cell->child[i];
                  if (bead->colour == 0) blue++;
                  if (bead->colour == 7) blac++;
                }
        if (blac) {
            cell->idata[2] = 0;	// being held in a bad contact
            cell->vdata[1] = cell->xyz; // record excited start position
        } else {
            if (blue>0) {		// running free with bad memories
                cell->idata[2]++;
                dd = cell->vdata[1] || cell->xyz;
                sumdd[cell->idata[2]] = dd;
                countMSD[cell->idata[2]]++;
            }
        }
        /*
if (cell->idata[2]) {
Pi(cell->id) Pi(cell->model) Pi(cell->idata[2]) Pr(dd) Pr(msd) Pv(cell->vdata[0]) NL
}
*/
    }
    DO(i,n) busyBead(cell->child[i]);	// implement behaviours
}

void setLinks ( Cell *bead ) {
    // set new links for the current <bead> based on <hit> found in bumper()
    int	it, me = bead->model;
    Data	*param = Data::model+me;
    int	level  = bead->level,
            nlinks = param->links[level],
            leader = bead->busy,	// WAS idata[0]
            hook, limit = 3;	// limit on the number of incoming links
    float	snap = param->sizes[level-1]; // break links over a bead diameter
    int	celhooks, hithooks;
    float unhook, collapse;
    Cell	*hit = bead->hit;	// the bead has bumpped another

    // ADD / CHANGE links between cells
    DO(i,nlinks)		// loop over links to check length
    { 
        Cell	*link = bead->link[i].to;
        float d;
        if (link==0) continue;
        d = bead->xyz | link->xyz;
        it = link->model;
        unhook = (Unhook*param->adh[abs(it - me)+1]);
        collapse = param->cil[abs(it - me)+1];
        if (Unhook>NOISE && randf()<unhook)
        {
            if (turns>0)
            {
                if (randf()<collapse)
                {
                    newLeader(bead->parent,turns);
                }
            }
            d = 999.0; // ie: break
        }
        if (d < snap) {
            continue;
        } else	// too far apart so break and decrement link counts
        { int	me1 = me+1, it1 = link->model+1;
            bead->nlinks--;		// uncount link
            bead->idata[it1]--;	// uncount incoming link to <bead>
            link->idata[me1]--;	// uncount incoming link to its <link>
            bead->link[i].to = 0;	// kill the link
            //Pt(kill) printf("%d.%d--%d.%d\n",bead->parent->id,bead->id,link->parent->id,link->id);
        }
    }
    if (hit==0) return; // no new interaction
    if (bead->parent == hit->parent) return; // an internal-hit
    celhooks = bead->idata[1]+ bead->idata[2]+ bead->idata[3];
    hithooks =  hit->idata[1]+  hit->idata[2]+  hit->idata[3];
    if (celhooks==limit || hithooks==limit) return;		// no room to hook
    DO(i,nlinks) { if (bead->link[i].to == hit) return; }	// check if already linked
    it = hit->model;
    // set <hook> according to sticking behaviour <stick> and bead types (<me> and <it>)
    hook = 0;
    switch (hooks) {	// 1=all, 2=RR+GG, 3=RR, 4=GG, 5=RG
    case 1 : hook = 1; break;
    case 2 : if (me==it) hook = 1; break;
    case 3 : if (me==it && me==1) hook = 1; break;
    case 4 : if (me==it && me==0) hook = 1; break;
    case 5 : if (me!=it) hook = 1; break;
    }
    // <hook>=1 means that the beads can stick to each other
    if (hook==0) return;
    // passed all checks so find a slot for the new link
    DO(i,nlinks)	// loop over links
    {
        Cell	*link = bead->link[i].to;
        int	me1 = me+1, it1 = it+1;
        if (link) continue;
        // found a free link slot
        link = bead->link[i].to = hit;	// new link made
        bead->nlinks++;			// count home out-links
        bead->idata[it1]++;	// count incoming link to this <bead>
        link->idata[me1]++;	// count incoming link to its <link>
        //Pt(make) printf("%d.%d--%d.%d\n",bead->parent->id,bead->id,link->parent->id,link->id);
        break;
    }
    bead->bump = 0; bead->hit = 0;
}

void memories ( Cell *bead ) {
    // Update the bead linkage memory in fdata[] from idata[] counts
    int	it, me = bead->model;
    Data	*param = Data::model+me;
    int	level  = bead->level,
            nlinks = param->links[level],
            leader = bead->busy,
            col;
    float	bias = 0.0;
    DO(i,NDATA) {
        bead->fdata[i] *= Decay; // degrade all linkage memories
    }
    DO(i,nlinks)	// loop over links and ++/-- memory in fdata[]
    { Cell	*link = bead->link[i].to;
        float f;
        if (link==0) continue;
        it = link->model;
        // increase linkage memory:  0 = leader, 1..3 = cell types
        if (link->busy) bead->fdata[0] += 1.0;	// bead linked to any leader
        if (bead->busy) link->fdata[0] += 1.0;	// link linked to any leader
        bead->fdata[it+1] += 1.0;
        link->fdata[me+1] += 1.0;
        f = bead->fdata[it+1];
        if (it==me) bias -= f; else bias += f;	// bad-good memories (diff-like)
    }
    // <bias> is the current total of bad-good memories
    bead->fdata[4] += bias;
    col = -1;	// default = same as cell
    if (bead->busy) { // colour leader yellow(5)/white(8) (excited parent)
        if (bead->parent->idata[2]) col = 8; else col = 5;
    } else {
        if (Flow && bead->fdata[0] > 1.0) {	// colour magenta(3)/cyan(4) = leader memory red/green(cyan)
            if (bead->model) col = 4; else col = 3;
        }
        if (bead->fdata[4] > 1.0) col = 0; 	// colour   blue(0) = bad memory
        if (bias > 1.0) col = 7;	// colour   black(7) = active bad contact
        //if (bias > 0.1) col = 0; 	// colour   blue(0) = bad memory
        //if (bias > 3.0) col = 7;	// colour   black(7) = very bad memory
    }
    // <colour> sets colour to be rendered in viewer()
    bead->colour = col;
    // 0=blue 1=green 2=red 3=cyan 4=magenta 5=yellow 6=grey 7=black 8=white
}

void busyBead ( Cell *bead ) {
    // animate the bead (cell sub-body)
    int	it, me = bead->model;
    Data	*param = Data::model+me;
    int	level  = bead->level,
            nlinks = param->links[level],
            excited = 0;
    bead->xyz.z *= 0.9;
    DO(i,nlinks)	// loop over links (currently just 1/bead)
    { Cell	*link = bead->link[i].to;
        if (link==0) continue;
        // found a link, so activate any associated behaviour
        it = link->model;
        if (bolts) {	// green runs away fast from red
            if ( me==0 && it==1) { // 0 = green
                bead->parent->fdata[0] = Fast;	// reset <speed> factor to panic value
                excited = 1;
            }
        }
        if (flees) {	// green moves away with greater persistence
           if ( me==0 && it==1) { // 0 = green
//            if ( me==0 ) { // 0 = green
                bead->parent->fdata[1] = Dive;	// reset <dodge> factor to panic value
                excited = 1;
            }
        }
    }
}

// penalties for interactions of type (big = don't select)
#define DIFF  3.0 //  20.0 //	 10.0	// different cell type
#define SAME  1.0 //   0.0 //	  0.0	// same cell type
#define LEAD 10.0 //  10.0 //    10.0   // memory of recent leader links
#define SEEN  5.0 //   5.0 //	  5.0 	// memory of recent alien links
#define STAY  0.0 //   2.0 //	  1.0	// bias to remain close to current
#define RAND  1.0 //   1.0 //	  1.0	// random noise (must not be zero)
#define EDGE 10.0 //  10.0 //	  1.0   // bias away from edge of the world

int newLeader ( Cell *cell, int avoid ) {
    // select a new leader bead: default = random
    // if <avoid> is set, then for cells of type <avoid>-1 (or all if > nmodels), select a leader
    // that is not linked to a cell of a different colour with preference given to unlinked cells
    // <avoid>=999 activates <EDGE> bias to avoid beads close to the edge of the world
    int	n = cell->kids,
            m = cell->model,
            links = Data::model[m].links[2],	// number of allotted bead links
            newboss, oldboss = cell->busy,	// leader number held in data[0]
            col;
    float	fear = 0.0;
    float	*bias = new float[n*2+1], f = twoPI/(float)n;
    float	*score = new float[n];
    int	*place = new int[n];
    bias[0] = 0.0;
    bias[n] = 1.0; // current bead position in ring
    DO(i,n) bias[n-i] = bias[n+i] = 1.0-0.5*(1.0+cos(f*(float)i)); // penalty to be far
    cell->child[oldboss]->busy =  0;
    cell->child[oldboss]->colour = -1;
    DO(i,n) // loop over beads to set basic behaviour (STAY near last, RANDom ,EDGE repulsion)
    { Cell *ci = cell->child[i];
        int	away = n+i-oldboss;		// (cyclic) displacement from leader
        score[i] = STAY * bias[away]	// bias to keep close to current position
                + RAND * randf();	// random addition
        if (avoid==999) { // at the edge of the world
            score[i] += EDGE * ci->xyz.len(); // add distance from centre
        }
    }
    if (avoid && links && (m+1==avoid || avoid>Data::nmodels)) { // add repeling behaviour (CI)
        DO(i,n) // loop over beads
        { Cell *ci = cell->child[i];
            float	same=0, diff=0;
    	    // other way to find linked cells
	    if (ci->nlinks)
	    {
		Cell	*link = ci->link[0].to;
		if (link->model==m) same++;
		else diff++;
	    }

            score[i] += SAME*same;		// links of like type
            score[i] += DIFF*diff;		// links of diff type

            score[(8+(i-1))%8] += 1.0/2.0*SAME*same;
            score[(8+(i-1))%8] += 1.0/2.0*DIFF*diff;
            score[(8+(i+1))%8] += 1.0/2.0*SAME*same;
            score[(8+(i+1))%8] += 1.0/2.0*DIFF*diff;
            fear += diff;
            if (memory) 			// memory of...
            { float lead = log(1.0+ci->fdata[0]),	// leader links (any sort)
                        seen = ci->fdata[4]; // [4] = bad-good memories, then take log difference
                if (seen > 0.0) seen = log(1.0+seen); else seen = -log(1.0-seen);
                score[i] += SEEN*seen;		// bad encounters
                if (Flow) score[i] += LEAD*lead; // and leader links
            }
        }
    }
    sort(score,place,n);
    newboss = place[n-1];
    cell->fdata[2] = score[place[0]];
    cell->fdata[3] = score[place[n-1]];
    cell->fdata[4] = cell->fdata[2]-cell->fdata[3];

    col = 5;	// default leader colour = yellow
    if (cell->idata[2]) col = 8;	// excited = white
    cell->child[newboss]->colour = col;
    cell->child[newboss]->busy = 1;	// mark leader bead
    cell->busy = newboss;	// record leader bead
    return newboss;
}
