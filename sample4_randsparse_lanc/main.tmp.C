//#include <alg/a2a/a2a.h>
//#include <alg/a2a/utils_main.h>
//#include <alg/a2a/grid_wrappers.h>
#include <cps.h>


USING_NAMESPACE_CPS

//A propagator G(x,x) = \sum_i V_i(x) W^\dag_i(x)
//Fortunately this is node local
//

CommonArg common_arg;
DoArg do_arg;
QPropWArg lqpropw_arg;
QPropWArg cqpropw_arg;
//LanczosArg lanczos_arg;
//MobiusArg mobius_arg;
MeasArg meas_arg;
//DoArgExt do_ext;
LancArg lanc_arg;
A2AArg a2a_arg; //note: src_width is the number of timeslices upon which the random source lives. To remove time dilution set src_width = Lt
CGcontrols cg_arg; //controls for the CG used to form the V vectors

inline int glb_sum(Rcomplex *c, const int n) {
  return glb_sum((Float*)c, 2*n);
}
inline int glb_sum_WMat(WilsonMatrix *W) {
  return glb_sum((Float*)W, 288);
}

Rcomplex get_F4(const WilsonMatrix &S1, const WilsonMatrix &S2,
                const WilsonMatrix &S3, const WilsonMatrix &S4)
{
  // Tr[ S1 S2 S3 S4 ]
  WilsonMatrix tmp1 = S1 * S2;
  WilsonMatrix tmp2 = S3 * S4;
  return Trace(tmp1,tmp2);
}

Rcomplex get_F3(const WilsonMatrix &S1, const WilsonMatrix &S2,
                const WilsonMatrix &S3)
{
  // Tr [ S1 S2 S3 ]
  WilsonMatrix tmp = S1 * S2;
  return Trace(tmp,S3);
}

//x is a *node local* lattice size 
template<typename A2Apolicies>
WilsonMatrix selfContractionBasic(int x[4], const A2AvectorV<A2Apolicies> &V, const A2AvectorW<A2Apolicies> &W){
  WilsonMatrix out;
  memset( (void*)&out, 0, 12*12*2*sizeof(double));

  int x3d = x[0] + GJP.XnodeSites()*( x[1] + GJP.YnodeSites()*x[2] );

  const int nmode = V.getNmodes();

  Rcomplex Vx[12], Wdagx[12];
  for(int mode=0;mode<nmode;mode++){
    for(int sc=0;sc<12;sc++){
      Vx[sc] = V.elem(mode, x3d, x[3], sc, 0);
      Wdagx[sc] = std::conj(W.elem(mode, x3d, x[3], sc, 0));
    }
    for(int s1=0;s1<4;s1++){
      for(int c1=0;c1<3;c1++){
	int sc1 = c1 + 3*s1;

	for(int s2=0;s2<4;s2++){
	  for(int c2=0;c2<3;c2++){
	    int sc2 = c2 + 3*s2;
	    
	    out(s1,c1,s2,c2) += Vx[sc1] * Wdagx[sc2];
	  }
	}
      }
    }
  }    
  return out;
}

//This should be tested by running with 0 low modes and using the high mode sources contained in the W field to generate normal volume stochastic propagators
//QPropWRandVolSrc uses the same 4D complex field for all spin and color
template<typename A2Apolicies>
void testRandVol(const A2AvectorW<A2Apolicies> &W, Lattice &lat, double mass){
//  assert(W.getNl() == 0); //makes sure 0 low modes
  assert(W.getNhits() == 1);
  typedef CPScomplex4D<cps::Rcomplex, FourDpolicy<DynamicFlavorPolicy>, StandardAllocPolicy> ComplexFieldType;
  const ComplexFieldType &rand_field = W.getWh(0);

  CgArg &cg = lqpropw_arg.cg;
  cg.mass = mass;
  cg.epsilon = 0;
  cg.max_num_iter = 10000;
  cg.stop_rsd = 1e-08;
  cg.true_rsd = 1e-08;
  cg.RitzMatOper = MATPCDAG_MATPC;

  lqpropw_arg.file = "";
  lqpropw_arg.flavor = 0;
  lqpropw_arg.gauge_fix_src = 0;
  lqpropw_arg.gauge_fix_snk = 0;
  lqpropw_arg.store_midprop = 0;
  lqpropw_arg.save_ls_prop = 0;
  lqpropw_arg.do_half_fermion = 0;
  lqpropw_arg.ensemble_label = "";
  lqpropw_arg.ensemble_id = "";
  lqpropw_arg.StartSrcSpin = 0;
  lqpropw_arg.EndSrcSpin = 3;
  lqpropw_arg.StartSrcColor = 0;
  lqpropw_arg.EndSrcColor = 2;

  CommonArg carg;
  QPropWRandArg rand_arg;
  
  QPropWRandVolSrc prop(lat, &lqpropw_arg, &rand_arg, &carg, (Complex const*)rand_field.ptr());
  
  if(!UniqueID()){
    int x[4] = {0,0,0,0};
    int site = x[0] + GJP.XnodeSites()*(x[1] + GJP.YnodeSites()*(x[2] + GJP.ZnodeSites()*x[3]));
    WilsonMatrix loop = prop[site];
    Complex wdag = std::conj(*rand_field.site_ptr(x));
    loop *= wdag;
    
    Rcomplex tr = loop.Trace();

    std::cout << "TEST: " << std::real(tr) << " " << std::imag(tr) << std::endl;
  }

}


//Defines a number of types used internally to the A2A library including the Grid fields and Dirac operator type
typedef A2ApoliciesDoubleAutoAlloc A2Apolicies;

int main (int argc,char **argv )
{
  const char *fname="main(int,char**)";
  Start(&argc, &argv);

  CommonArg carg;

  if(!do_arg.Decode("do_arg.vml","do_arg")){
    do_arg.Encode("do_arg.dat","do_arg");
    VRB.Result("","main","Can't open do_arg.vml!\n");exit(1);
  }

  if(!lanc_arg.Decode("lanc_arg.vml","lanc_arg")){
    lanc_arg.Encode("lanc_arg.templ","lanc_arg");
    VRB.Result("","main","Can't open lanc_arg.vml!\n");exit(1);
  }
  if(!a2a_arg.Decode("a2a_arg.vml","a2a_arg")){
    a2a_arg.Encode("a2a_arg.templ","a2a_arg");
    VRB.Result("","main","Can't open a2a_arg.vml!\n");exit(1);
  }

  if(!cg_arg.Decode("cg_arg.vml","cg_arg")){
    cg_arg.Encode("cg_arg.templ","cg_arg");
    VRB.Result("","main","Can't open cg_arg.vml!\n");exit(1);
  }

  if(!meas_arg.Decode("meas_arg.vml","meas_arg")){
    meas_arg.Encode("meas_arg.templ","meas_arg");
    VRB.Result("","main","Can't open meas_arg.vml!\n");exit(1);
  }

  if(!lqpropw_arg.Decode("lqpropw_arg.vml","lqpropw_arg")){
    lqpropw_arg.Encode("lqpropw_arg.templ","lqpropw_arg");
    VRB.Result("","main","Can't open lqpropw_arg.vml!\n");exit(1);
  }
  std::cout << "mass of light quarks: " << lqpropw_arg.cg.mass << std::endl;

  if(!cqpropw_arg.Decode("cqpropw_arg.vml","cqpropw_arg")){
    cqpropw_arg.Encode("cqpropw_arg.templ","cqpropw_arg");
    VRB.Result("","main","Can't open cqpropw_arg.vml!\n");exit(1);
  }
  std::cout << "mass of charm quark: " << cqpropw_arg.cg.mass << std::endl;

  //Read command line args                                                   

  int N_evec= 0;
  N_evec= lanc_arg.N_true_get;


  //begin to work
  chdir(meas_arg.WorkDirectory);

  //main loop for each configuration
  /*
  int low = meas_arg.TrajStart;
  int high = meas_arg.TrajLessThanLimit;
  int dt = meas_arg.TrajIncrement;
  */
  int low = atoi(argv[2]);
  int high = atoi(argv[3]);
  int dt = atoi(argv[4]);


  const int nCG_psrc = atoi(argv[5]);
  const int nsrc_1CG = atoi(argv[6]);
  const int nsrc = nCG_psrc * nsrc_1CG;
  const int nr = atoi(argv[7]);
  std::cout << "# of CG invs for point sources: " << nCG_psrc << std::endl;
  std::cout << "Npsrc per CG: " << nsrc_1CG << std::endl;
  std::cout << "Npsrc: " << nsrc << std::endl;


  //Setup CPS
  //do_arg.start_seed_value += low+high + nr*840923 + 137*nsrc;
  do_arg.start_seed_value += low*137;

  GJP.Initialize(do_arg);
//  GJP.InitializeExt(do_ext);

  std::cout << "Seed: " << GJP.StartSeedValue() << std::endl;
  LRG.Initialize();

  //int nthreads = GJP.Nthreads();
  int nthreads = omp_get_max_threads();
#ifdef USE_GRID
  std::cout << "Number of threads " << nthreads << std::endl;
  Grid::GridThread::SetThreads(nthreads);
#endif


  std::cout << "start: " << low << ",  max: " << high << ",  step: " << dt << std::endl;

  QPropWRandArg rarg;
  rarg.rng = UONE;
  rarg.seed = 1234*low;// !
  for(meas_arg.TrajCur = low; meas_arg.TrajCur < high; meas_arg.TrajCur += dt) {
    rarg.seed += meas_arg.TrajCur;
    std::cout << "configuration " << meas_arg.TrajCur << " start :" << std::endl;
    double cost = -dclock();

    ReadGaugeField(meas_arg);
  
    //Initialize FGrid
    FgridParams fgp;

    //Mobius parameters
    const double mob_b = 1.0;
    const double mob_c = 0.0;   //b-c = 1
    fgp.mobius_scale = mob_b + mob_c; //b+c
    const double M5 = do_arg.dwf_height;
    printf("Grid b=%g c=%g b+c=%g\n",mob_b,mob_c,mob_b+mob_c);

    typedef A2Apolicies::FgridGFclass LatticeType;
    typedef typename A2Apolicies::GridFermionField GridFermionField;
    typedef typename A2Apolicies::GridFermionFieldF GridFermionFieldF;
    typedef typename A2Apolicies::FgridFclass FgridFclass;
    typedef typename A2Apolicies::GridDirac GridDirac;
    LatticeType lattice(fgp);
    lattice.ImportGauge();

    //Do the Lanczos
    std::cout << "GridLanczosWrapper<A2Apolicies> lanczos" << endl;
    GridLanczosWrapper<A2Apolicies> lanczos;
    if(N_evec>0){
      std::cout << "lanczos.compute(lanc_arg, lattice)" << endl;
      lanczos.compute(lanc_arg, lattice);
    }

    //Typically we convert the evecs to single precision to save memory
    //(Note split Grid not yet supported if we don't have single-prec evecs)
    lanczos.toSingle(lattice);
  
    A2AvectorW<A2Apolicies> W(a2a_arg);
    A2AvectorV<A2Apolicies> V(a2a_arg);

    if(N_evec>0)
      W.computeVW(V, lattice, lanczos.evec_f, lanczos.eval, lanc_arg.mass, cg_arg);

    EvecInterfaceGridSinglePrec<A2Apolicies> ev(lanczos.evec_f,lanczos.eval,lattice,lanc_arg.mass);

    FgridFclass &latg = dynamic_cast<FgridFclass&>(lattice);
#if 0
    //Grids and gauge field
    Grid::GridCartesian *UGrid = latg.getUGrid();
    Grid::GridRedBlackCartesian *UrbGrid = latg.getUrbGrid();
    Grid::GridCartesian *FGrid = latg.getFGrid();
    Grid::GridRedBlackCartesian *FrbGrid = latg.getFrbGrid();
    Grid::LatticeGaugeFieldD *Umu = latg.getUmu();
    std::cout<<"UGrid= "<<UGrid<<" UrbGrid= "<<UrbGrid<<" FGrid= "<<FGrid<<" FrbGrid= "<<FrbGrid<<std::endl;
#endif

    //const int gparity = GJP.Gparity();

    //Setup Grid Dirac operator
    typename GridDirac::ImplParams params;
    latg.SetParams(params);


    if(!UniqueID()){
      int x[4] = {0,0,0,0};
      WilsonMatrix loop = selfContractionBasic(x,V,W);
    
      Rcomplex tr = loop.Trace();

      std::cout << "Loop trace, site 0: " << std::real(tr) << " " << std::imag(tr) << std::endl;
    }


    const char *evec_name = "light_evec";

    lqpropw_arg.cg.fname_eigen = (char *) evec_name;

    lqpropw_arg.cg.neig=N_evec;

    if(N_evec>0)
    {
      EigenCacheGrid<GridFermionFieldF>  *ecache = new EigenCacheGrid <GridFermionFieldF> (evec_name);
//    const int n_fields = GJP.SnodeSites ();
//    const size_t f_size_per_site = lattice.FsiteSize () / n_fields / 2;     // checkerboarding
//    size_t evec_size = (size_t) (GJP.VolNodeSites () / 2) * lattice.FsiteSize ();
//    assert(evec_size != lattice.half_size)
//    size_t fsize = evec_size;
//    int data_size = sizeof (Float);
//    if (lanczos_arg.precision == PREC_SINGLE)
//      data_size = sizeof (float);
//    ecache->alloc (N_evec,evec_size, data_size);
      ecache->load(lanczos.eval, lanczos.evec_f);
//    ecache->read_compressed ((char*)evec_dir);
      EigenCacheList.push_back (ecache);
      lqpropw_arg.cg.Inverter = CG_LOWMODE_DEFL;
    } else {
      lqpropw_arg.cg.Inverter = CG;
    }

    std::cout << "Seed: " << GJP.StartSeedValue() << std::endl;
    LRG.Reinitialize();

    std::cout << "Nnoise: " << nr << std::endl;
    std::cout << "Seed: " << do_arg.start_seed_value << std::endl;
    std::cout << "Start inversion with 1st noise (light)" << std::endl;
    QPropWRandVolSrc rnd0(lattice, &lqpropw_arg, &rarg, &common_arg);
    std::vector<QPropWRandVolSrc> lrprop(nr,rnd0);
    for( int ir = 1; ir < nr; ++ir ) {
      std::cout << "Start inversion with " << (ir+1) << "th noise (light)" << std::endl;
      QPropWRandVolSrc rnd(lattice, &lqpropw_arg, &rarg, &common_arg);
      lrprop[ir] = rnd;
    }

    rarg.seed += 4791 * nr;
    std::cout << "Seed: " << do_arg.start_seed_value << std::endl;
    std::cout << "Nnoise: " << nr << std::endl;
    std::cout << "Start inversion with 1st noise (charm)" << std::endl;
    QPropWRandVolSrc rnd1(lattice, &cqpropw_arg, &rarg, &common_arg);
    std::vector<QPropWRandVolSrc> crprop(nr,rnd1);
    for( int ir = 1; ir < nr; ++ir ) {
      std::cout << "Start inversion with " << (ir+1) << "th noise (charm)" << std::endl;
      QPropWRandVolSrc rnd(lattice, &cqpropw_arg, &rarg, &common_arg);
      crprop[ir] = rnd;
    }

    int mu,nu;
    const int glb_x = GJP.XnodeSites() * GJP.Xnodes();
    const int glb_y = GJP.YnodeSites() * GJP.Ynodes();
    const int glb_z = GJP.ZnodeSites() * GJP.Znodes();
    const int glb_t = GJP.TnodeSites() * GJP.Tnodes();
    const long long int glb_v
      = GJP.XnodeSites() * GJP.Xnodes()
      * GJP.YnodeSites() * GJP.Ynodes()
      * GJP.ZnodeSites() * GJP.Znodes()
      * GJP.TnodeSites() * GJP.Tnodes();

    // __VECTOR_FORMATS__
    std::vector<Rcomplex> corr_pp(glb_t, Rcomplex(0.,0.));
    std::vector<Rcomplex> corr_disc(glb_t, Rcomplex(0.,0.));

    for (int isrc = 0; isrc < nCG_psrc; ++isrc ) {
      int src_x,src_y,src_z,src_t;
      if ( nsrc_1CG == 1 ) {
	if ( nsrc == 64 ) {
	  int id = isrc % 2;
	  int tmp_isrc = int(isrc / 2);
	  int ix = tmp_isrc % 2;
	  tmp_isrc = int(tmp_isrc / 2);
	  int iy = tmp_isrc % 2;
	  tmp_isrc = int(tmp_isrc / 2);
	  int iz = tmp_isrc % 2;
	  int it = int(tmp_isrc / 2);
	  src_x = (2*ix + id) * glb_x/4;
	  src_y = (2*iy + id) * glb_y/4;
	  src_z = (2*iz + id) * glb_z/4;
	  src_t = (2*it + id) * glb_t/8;
	} else if ( nsrc == 16 ) {
	  int isrc_s = isrc % 4;
	  int isrc_t = int (isrc / 4);
	  int ix = int(isrc_s / 2);
	  int iy = isrc_s % 2;
	  int iz = 0;
	  if ( ix != iy ) {
	    ++ix; ++iy; ++iz;
	    ix = ix % 2; iy = iy % 2; iz = iz % 2;
	  }
	  ix = 2 * ix + isrc_t; iy = 2 * iy + isrc_t; iz = 2 * iz + isrc_t;
	  ix = ix % 4; iy = iy % 4; iz = iz % 4;
	  src_x = ix * glb_x/4;
	  src_y = iy * glb_y/4;
	  src_z = iz * glb_z/4;
	  src_t = isrc_t*glb_t/4;
	} else if ( nsrc == 4 ) {
	  int iz = isrc % 2;
	  int ix = iz;
	  int iy = iz;
	  src_x = ix * glb_x/2;
	  src_y = iy * glb_y/2;
	  src_z = iz * glb_z/2;
	  src_t = isrc*glb_t/4;
	} else if ( nsrc == 2 || nsrc == 1 ) {
	  src_x = isrc * glb_x/2;
	  src_y = isrc * glb_y/2;
	  src_z = isrc * glb_z/2;
	  src_t = isrc * glb_t/2;
	} else {
          std::cout << "Source Points with nsrc = "
		    << nsrc << " = " << nsrc_1CG << " x " << nCG_psrc
		    << " not registered\n" << std::endl;
	}
      } else if ( nsrc_1CG == 4 ) {
        if ( nCG_psrc == 16 ) {
          int id = int(isrc/8);
          int ixyz = isrc%8;
          int iz = int(ixyz/4);
          int ixy = ixyz%4;
          int iy = int(ixy/2);
          int ix = ixy%2;
          src_x = ( 2*ix + id ) * glb_x/4;
          src_y = ( 2*iy + id ) * glb_y/4;
          src_z = ( 2*iz + id ) * glb_z/4;
          src_t = id * glb_t/8;
        }
      }
      lqpropw_arg.x = src_x;
      lqpropw_arg.y = src_y;
      lqpropw_arg.z = src_z;
      lqpropw_arg.t = src_t;
      cqpropw_arg.x = src_x;
      cqpropw_arg.y = src_y;
      cqpropw_arg.z = src_z;
      cqpropw_arg.t = src_t;

      std::cout << "Start lquark inversions with " << (isrc+1) << "th source point: "
		<< lqpropw_arg.x << " " << lqpropw_arg.y << " "
		<< lqpropw_arg.z << " " << lqpropw_arg.t << std::endl;
      std::cout << "1st light source" << std::endl;
      QPropWRandArg sparse_arg;
      sparse_arg.rng=ZTHREE;
      sparse_arg.sep=nsrc_1CG;
      QPropWRandSparse lqpropw1(lattice, &lqpropw_arg, &sparse_arg, &common_arg);
      std::cout << "2nd light source" << std::endl;
      QPropWRandSparse lqpropw2(lattice, &lqpropw_arg, &sparse_arg, &common_arg);

      std::cout << "Start cquark inversion with " << (isrc+1) << "th source point: "
		<< cqpropw_arg.x << " " << cqpropw_arg.y << " "
		<< cqpropw_arg.z << " " << cqpropw_arg.t << std::endl;
      QPropWRandSparse cqpropw(lattice, &cqpropw_arg, &sparse_arg, &common_arg);
      //QPropWPointSrc cqpropw(lattice, &cqpropw_arg, &common_arg);

      for ( int iy0 = 0; iy0 < nsrc_1CG; ++iy0 ) {
	int x0, y0, z0, t0;
	lqpropw1.SourcePointsRel(iy0,nsrc_1CG,x0,y0,z0,t0);
	x0 = ( x0 + src_x ) % glb_x;
	y0 = ( y0 + src_y ) % glb_y;
	z0 = ( z0 + src_z ) % glb_z;
	t0 = ( t0 + src_t ) % glb_t;
	std::cout << (iy0+1) << "th source point: "
		  << x0 << " " << y0 << " " << z0 << " " << t0 << std::endl;

	WilsonMatrix lpropyy1 = lqpropw1[0];
	WilsonMatrix lpropyy2 = lqpropw2[0];
	WilsonMatrix cpropyy = cqpropw[0];// ! Will be changed so that each node has the same quark-loop propagator at the source point
	lpropyy1 *= 0.0;
	lpropyy2 *= 0.0;
	cpropyy *= 0.0;

	Site s0(0);
	Rcomplex xi_l1 = 0.,xi_l2 = 0.,xi_c=0.;
	while(s0.End()){
	  int xsq = ( s0.physX() - x0 ) * ( s0.physX() - x0 )
	    +       ( s0.physY() - y0 ) * ( s0.physY() - y0 )
	    +       ( s0.physZ() - z0 ) * ( s0.physZ() - z0 )
	    +       ( s0.physT() - t0 ) * ( s0.physT() - t0 );
	  if ( xsq > 0 ) {
	    s0.nextSite();
	    continue;
	  }
	  xi_l1 = Rcomplex(lqpropw1.rsrc[2*s0.Index()],
			   lqpropw1.rsrc[2*s0.Index()+1]);
	  xi_l2 = Rcomplex(lqpropw2.rsrc[2*s0.Index()],
			   lqpropw2.rsrc[2*s0.Index()+1]);
	  xi_c  = Rcomplex(cqpropw.rsrc[2*s0.Index()],
			   cqpropw.rsrc[2*s0.Index()+1]);
	  lpropyy1 = lqpropw1[s0.Index()];
	  lpropyy2 = lqpropw2[s0.Index()];
	  cpropyy = cqpropw[s0.Index()];
	  s0.nextSite();
	}// s0
	glb_sum(&xi_l1,1);
	glb_sum(&xi_l2,1);
	glb_sum(&xi_c,1);
	glb_sum_WMat(&lpropyy1);
	glb_sum_WMat(&lpropyy2);
	glb_sum_WMat(&cpropyy);

	// __TRACES_AT_Y__
	Rcomplex F1yg5l(0.,0.);
	{
	  WilsonMatrix tmp = lpropyy1 * conj(xi_l1);
	  tmp = tmp.gl(-5);
	  F1yg5l = Trace(tmp);
	  tmp = lpropyy2 * conj(xi_l2);
	  tmp = tmp.gl(-5);
	  F1yg5l += Trace(tmp);
	  F1yg5l *= 0.5;
	}

	Site s(0);
	while(s.End()){ //Logically this should be !s.End() but whoever wrote the code for Site was clearly lacking the appropriate amount of caffeine
	  int rel_x, rel_y, rel_z, rel_t;
	  rel_x = (s.physX() + glb_x - x0) % glb_x;
	  rel_y = (s.physY() + glb_y - y0) % glb_y;
	  rel_z = (s.physZ() + glb_z - z0) % glb_z;
	  rel_t = (s.physT() + glb_t - t0) % glb_t;
	  long long int id
	    = rel_x+glb_x*(rel_y+glb_y*(rel_z+glb_z*rel_t));

	  WilsonMatrix lpropxy1 = lqpropw1[s.Index()];
	  WilsonMatrix lpropxy2 = lqpropw2[s.Index()];
	  WilsonMatrix lpropyx1 = lpropxy1;
	  lpropyx1.hconj();
	  WilsonMatrix lpropyx2 = lpropxy2;
	  lpropyx2.hconj();

	  Rcomplex r = Trace(lpropxy1,lpropyx2);  //this is a simple point-source pseudoscalar correlation function
	  corr_pp[rel_t] += r * conj(xi_l1) * xi_l2;

	  // gamma_5 sandwich
	  lpropyx1.gr(-5); lpropyx1.gl(-5);
	  lpropyx2.gr(-5); lpropyx2.gl(-5);

	  WilsonMatrix lpropxx = lpropxy1;
	  lpropxx *= 0;
	  WilsonMatrix cpropxx = lpropxx;
	  for(int ir = 0; ir < nr; ++ir) {
	    lpropxx +=
	      lrprop[ir][s.Index()] * conj(lrprop[ir].rand_src(s.Index()));
	    cpropxx +=
	      crprop[ir][s.Index()] * conj(crprop[ir].rand_src(s.Index()));
	    //= rnd[s.Index()] * conj(rnd.rand_src(s.Index()));
	  }
	  lpropxx *= double(1./nr);
	  cpropxx *= double(1./nr);

	  WilsonMatrix propx1; // will be Γx1 S(x-[xy])
	  WilsonMatrix propx2; // will be Γx2 S(x-[xy])
	  WilsonMatrix propy1; // will be Γy1 S(y-[xy])
	  WilsonMatrix propy2; // will be Γy2 S(y-[xy])

	  // __TRACES_AT_X_ANDOR_Y__
	  Rcomplex F1xg5l(0.,0.);
	  propx1 = lpropxx;
	  propx1 = propx1.gl(-5);
	  F1xg5l = Trace(propx1);

	  for ( mu = 0; mu <= 3; ++mu ) {
	    for ( nu = 0; nu <= 3; ++nu ) {
	      int munu = 4*mu + nu;
	      // __ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_MUNU__
	    }
	    // __ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_MU__
	  }
	  // __ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_SCL__
	  corr_disc[rel_t] += F1xg5l * F1yg5l;
	  s.nextSite();
	}// loop over sites
      }// loop over source points iy0 < nsrc_1CG
    }// loop over isrc < nCG_psrc
    // __SUM_OVER_ALL_NODES__
    glb_sum(corr_pp.data(), corr_pp.size()); 
    glb_sum(corr_disc.data(), corr_pp.size());

    std::fstream file;
    Site s(0);
    int blk_x = s.physX() / GJP.XnodeSites();
    int blk_y = s.physY() / GJP.YnodeSites();
    int blk_z = s.physZ() / GJP.ZnodeSites();
    int blk_t = s.physT() / GJP.TnodeSites();
    int nodeid = blk_x + GJP.Xnodes() * ( blk_y + GJP.Ynodes() * ( blk_z + GJP.Znodes() * blk_t ));
    int nnode = GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
    // __OUTPUT_ON_BINARY_FILES__

    std::cout << "Pseudoscalar Corr:\n";
    for(int t=0;t<corr_pp.size();t++) std::cout << t << " " << corr_pp[t]/nsrc << " " << corr_disc[t]/nsrc << std::endl;

    cost += dclock();
    if (s.physX()==0 && s.physY()==0 && s.physZ()==0 && s.physT()==0 )
      std::printf("configuration %d end, time spend = %f sec. \n",meas_arg.TrajCur, cost);
  }
  std::cout << "The End!" << std::endl;

  End();
  return 0;
}



