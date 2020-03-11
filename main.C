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

#if 0
void ReadGaugeField(const MeasArg &meas_arg)
{
  const char *cname = "main";
  const char *fname = "ReadGaugeField";

  GnoneFnone lat;
  char lat_file[256];
  sprintf(lat_file, "%s.%d", meas_arg.GaugeStem, meas_arg.TrajCur);
  VRB.Result(cname, fname, lat_file);
  QioArg rd_arg(lat_file, 0.001);//chkprec=0.001
  rd_arg.ConcurIONumber=meas_arg.IOconcurrency;
  ReadLatticeParallel rl;
  rl.read(lat, rd_arg);
  if(!rl.good()) 
    ERR.General(cname, fname, "Failed read lattice %s", lat_file);	
}
#endif

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


  const int nsrc = atoi(argv[5]);
  const int nr = atoi(argv[6]);
  do_arg.start_seed_value += low+high + nr*840923 + 137*nsrc;


  //Setup CPS
  GJP.Initialize(do_arg);
//  GJP.InitializeExt(do_ext);
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
    std::vector<Rcomplex> CON_F4gmuLlgnuLlgmuLlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4gmuLlgLlgmuLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4gLlgnuLlgRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4gmuLlgnuLlgmuLlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4gmuLlgnuLlgmuRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4gmuLlgRlgmuRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4gLlgnuLlgRlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4gmuLlgnuLlgmuRlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4gLlgLlgRlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F2gmuLlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F2gmuLlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F2gmuLlgnuRlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F2gmuLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F2gRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F2gmuLlgnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgRl_F2gRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F2gmuLlgRlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F2gRlgnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuLlgnuLlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuLlgnuLcgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLcgmuLlgnuLlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuLlgRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLlgRlgnuRlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuLlgnuRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuRlgmuLlgnuLlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuLlgnuLlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuRlgnuLlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuLlgRcgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLcgRlgnuRlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLlgRlgnuRcgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLcgmuLlgRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLlgRlgLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuRlgmuLlgLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLlgRlgnuRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuRlgLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLlgRlgnuLlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLlgRlgLcgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLcgRlgLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuRlgnuLcgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLcgmuLlgnuLlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuRlgnuRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuRlgmuLlgnuLlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuRlgnuLlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLlgmuRlgLcgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLcgRlgnuLlgnuRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuRlgmuLlgnuLcgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuLcgmuLlgnuRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuRlgmuLlgnuRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgmuRlgmuLlgLcgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F4pgLcgRlgnuRlgnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuLlgnuLl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuLlgmuLl_F1xgmuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLcgmuLlgnuLl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLcgnuLlgmuLl_F1xgmuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLlgRlgnuRlC_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLlgRlgmuRlC_F1xgmuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuRlgmuLlgnuLl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuRlgnuLlgmuLl_F1xgmuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuRlgnuLl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuRlgmuLl_F1xgmuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLcgRlgnuRlC_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLcgRlgmuRlC_F1xgmuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuLlgnuLl_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuLlgmuLl_F1xgmuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuLlgnuLl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuLlgmuLl_F1xgmuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuLlgLl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuLlgLl_F1xgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuLlgRl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuLlgRl_F1xgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuLlgnuLl_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuLlgmuLl_F1xgmuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLcgmuLlgnuLl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLcgnuLlgmuLl_F1xgmuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLlgRlgnuRlC_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLlgRlgmuRlC_F1xgmuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuRlgmuLlgnuLl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuRlgnuLlgmuLl_F1xgmuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuRlgnuLl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuRlgmuLl_F1xgmuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLcgRlgnuRlC_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLcgRlgmuRlC_F1xgmuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLlgRlgnuRlC_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLlgRlgmuRlC_F1xgmuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLlgRlgLlC_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLlgRlgLlC_F1xgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLlgRlgRlC_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLlgRlgRlC_F1xgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLlgRlgnuRlC_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLlgRlgmuRlC_F1xgmuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuRlgnuLl_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuRlgmuLl_F1xgmuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuRlgLl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuRlgLl_F1xgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuRlgRl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuRlgRl_F1xgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuRlgnuLl_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuRlgmuLl_F1xgmuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLcgmuLlgRl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLcgnuLlgRl_F1xgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuRlgmuLlgRl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuRlgnuLlgRl_F1xgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLcgRlgRlC_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLcgRlgRlC_F1xgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuRlgmuLlgnuLl_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuRlgnuLlgmuLl_F1xgmuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuRlgmuLlgLl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuRlgnuLlgLl_F1xgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuRlgmuLlgnuLl_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuRlgnuLlgmuLl_F1xgmuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLcgmuLlgLl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLcgnuLlgLl_F1xgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLcgRlgLlC_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLcgRlgLlC_F1xgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F1xgmuLl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F1xgRl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgRl_F1xgmuLl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgnuLl_F1xgLl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F1xgmuLlC_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F1xgRl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgRl_F1xgmuLlC_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgnuLl_F1xgLl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgRl_F1xgmuLc_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgnuLl_F1xgLl_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgRl_F1xgRl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgLl_F1xgLl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgRl_F1xgLl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgRl_F1xgmuLcC_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgnuLl_F1xgLl_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F1xgmuLc_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F1xgRl_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgLl_F1xgRl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F1xgmuLcC_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F1xgRl_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuLlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuLlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLlgRlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLlgRlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLlgRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLlgRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuRlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuRlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLlgmuRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLlgnuRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuRlgmuLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuRlgnuLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuRlgmuLlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuRlgnuLlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLcgmuLlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLcgnuLlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLcgRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLcgRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgmuLcgmuLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygnuLcgnuLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3xgLcgRlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F3ygLcgRlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgRl_F1xgmuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgnuLl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F1xgmuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F1ygnuLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgRl_F1xgmuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgnuLl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F1xgmuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F1ygnuLlC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgRl_F1xgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgRl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgLl_F1xgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgRl_F1ygLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgRl_F1xgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgLl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgLl_F1xgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgLl_F1ygRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F1xgmuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgLl_F1xgmuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgnuLl_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgRl_F1xgmuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgnuLl_F1ygnuLc(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gmuLlgRl_F1xgmuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgnuLl_F1ygnuLcC(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gLlgLl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> CON_F2gRlgRl(glb_v, Rcomplex(0.,0.));
    std::vector<Rcomplex> corr_pp(glb_t, Rcomplex(0.,0.));
    std::vector<Rcomplex> corr_disc(glb_t, Rcomplex(0.,0.));

    std::cout << "Nsrc: " << nsrc << std::endl;
    for (int isrc = 0; isrc < nsrc; ++isrc ) {
      int src_x,src_y,src_z,src_t;
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
	std::cout << "Source Points with nsrc = " << nsrc <<
	  " not registered\n" << std::endl;
      }
      lqpropw_arg.x = src_x;
      lqpropw_arg.y = src_y;
      lqpropw_arg.z = src_z;
      lqpropw_arg.t = src_t;
      cqpropw_arg.x = src_x;
      cqpropw_arg.y = src_y;
      cqpropw_arg.z = src_z;
      cqpropw_arg.t = src_t;

      std::cout << "Start lquark inversion with " << isrc << "th source point: "
		<< lqpropw_arg.x << " " << lqpropw_arg.y << " "
		<< lqpropw_arg.z << " " << lqpropw_arg.t << std::endl;
      QPropWPointSrc lqpropw(lattice, &lqpropw_arg, &common_arg);

      std::cout << "Start cquark inversion with " << isrc << "th source point: "
		<< cqpropw_arg.x << " " << cqpropw_arg.y << " "
		<< cqpropw_arg.z << " " << cqpropw_arg.t << std::endl;
      QPropWPointSrc cqpropw(lattice, &cqpropw_arg, &common_arg);


      WilsonMatrix lpropyy = lqpropw[0];// ! Will be changed so that each node has the same quark-loop propagator at the source point
      WilsonMatrix cpropyy = cqpropw[0];// ! Will be changed so that each node has the same quark-loop propagator at the source point

      Site s0(0);
      int rel_x, rel_y, rel_z, rel_t;
      lpropyy *= 0.0;
      cpropyy *= 0.0;

      while(s0.End()){
        int x2 =
            ( s0.physX() - src_x ) * ( s0.physX() - src_x )
          + ( s0.physY() - src_y ) * ( s0.physY() - src_y )
          + ( s0.physZ() - src_z ) * ( s0.physZ() - src_z )
          + ( s0.physT() - src_t ) * ( s0.physT() - src_t );
        if ( x2 > 0 ) {
          s0.nextSite();
          continue;
        }
        lpropyy = lqpropw[s0.Index()];
        cpropyy = cqpropw[s0.Index()];
        printf("source loop taken from: %d %d %d %d\n",s0.physX(),s0.physY(),s0.physZ(),s0.physT());
        s0.nextSite();
      }
      glb_sum_WMat(&lpropyy);
      glb_sum_WMat(&cpropyy);

      // __TRACES_AT_Y__
      std::vector<Rcomplex> F1ygnuLl(4, Rcomplex(0.,0.));
      for ( nu = 0; nu <= 3; ++nu ) {
        WilsonMatrix tmp = lpropyy;
        tmp = tmp.glL(nu);// γν(1-γ5) S(y-y)
        F1ygnuLl[nu] = Trace(tmp);
      }
      std::vector<Rcomplex> F1ygnuLl(4, Rcomplex(0.,0.));
      for ( nu = 0; nu <= 3; ++nu ) {
        WilsonMatrix tmp = lpropyy;
        tmp = tmp.glL(nu);// γν(1-γ5) S(y-y)
        F1ygnuLl[nu] = Trace(tmp);
      }
      std::vector<Rcomplex> F1ygnuLc(4, Rcomplex(0.,0.));
      for ( nu = 0; nu <= 3; ++nu ) {
        WilsonMatrix tmp = cpropyy;
        tmp = tmp.glL(nu);// γν(1-γ5) S(y-y)
        F1ygnuLc[nu] = Trace(tmp);
      }
      Rcomplex F1ygRl(0.,0.);
      {
        WilsonMatrix tmp = lpropyy;
        tmp = tmp.glPR();// 0.5 (1+γ5) S(y-y)
        tmp *= 2.0;
        F1ygRl = Trace(tmp);
      }
      Rcomplex F1ygRl(0.,0.);
      {
        WilsonMatrix tmp = lpropyy;
        tmp = tmp.glPR();// 0.5 (1+γ5) S(y-y)
        tmp *= 2.0;
        F1ygRl = Trace(tmp);
      }
      Rcomplex F1ygLl(0.,0.);
      {
        WilsonMatrix tmp = lpropyy;
        tmp = tmp.glPL();// 0.5 (1-γ5) S(y-y)
        tmp *= 2.0;
        F1ygLl = Trace(tmp);
      }
      Rcomplex F1ygLl(0.,0.);
      {
        WilsonMatrix tmp = lpropyy;
        tmp = tmp.glPL();// 0.5 (1-γ5) S(y-y)
        tmp *= 2.0;
        F1ygLl = Trace(tmp);
      }
      Rcomplex F1yg5l(0.,0.);
      {
	WilsonMatrix tmp = lpropyy;
	tmp = tmp.gl(-5);
	F1yg5l = Trace(tmp);
      }

      Site s(0);
      while(s.End()){ //Logically this should be !s.End() but whoever wrote the code for Site was clearly lacking the appropriate amount of caffeine
	rel_x = (s.physX() + glb_x - src_x) % glb_x;
	rel_y = (s.physY() + glb_y - src_y) % glb_y;
	rel_z = (s.physZ() + glb_z - src_z) % glb_z;
	rel_t = (s.physT() + glb_t - src_t) % glb_t;
	long long int id
	  = rel_x+glb_x*(rel_y+glb_y*(rel_z+glb_z*rel_t));

	WilsonMatrix lpropxy = lqpropw[s.Index()];
	WilsonMatrix lpropyx = lpropxy;
	lpropyx.hconj();

	Rcomplex r = Trace(lpropxy,lpropyx);  //this is a simple point-source pseudoscalar correlation function
	corr_pp[rel_t] += r;

	lpropyx.gr(-5); //right-multiply by gamma^5
	lpropyx.gl(-5); //left-multiply by gamma^5

	WilsonMatrix cpropxy = cqpropw[s.Index()];
	WilsonMatrix cpropyx = cpropxy;
	cpropyx.hconj();
	cpropyx.gr(-5); //right-multiply by gamma^5
	cpropyx.gl(-5); //left-multiply by gamma^5

	WilsonMatrix lpropxx = lpropxy;
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
        Rcomplex F4gmuLlgnuLlgmuLlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4gmuLlgnuLlgmuLlgnuLl += get_F4(propx1,propy1,propx2,propy2);
          }// nu
        }// mu
        Rcomplex F4gmuLlgLlgmuLlgRl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F4gmuLlgLlgmuLlgRl += get_F4(propx1,propy1,propx2,propy2);
        }// mu
        Rcomplex F4gLlgnuLlgRlgnuLl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyx;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F4gLlgnuLlgRlgnuLl += get_F4(propx1,propy1,propx2,propy2);
        }// nu
        Rcomplex F4gmuLlgnuLlgmuLlgnuRl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glR(nu);
            F4gmuLlgnuLlgmuLlgnuRl += get_F4(propx1,propy1,propx2,propy2);
          }// nu
        }// mu
        Rcomplex F4gmuLlgnuLlgmuRlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4gmuLlgnuLlgmuRlgnuLl += get_F4(propx1,propy1,propx2,propy2);
          }// nu
        }// mu
        Rcomplex F4gmuLlgRlgmuRlgLl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPL();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          F4gmuLlgRlgmuRlgLl += get_F4(propx1,propy1,propx2,propy2);
        }// mu
        Rcomplex F4gLlgnuLlgRlgnuRl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyx;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          F4gLlgnuLlgRlgnuRl += get_F4(propx1,propy1,propx2,propy2);
        }// nu
        Rcomplex F4gmuLlgnuLlgmuRlgnuRl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glR(nu);
            F4gmuLlgnuLlgmuRlgnuRl += get_F4(propx1,propy1,propx2,propy2);
          }// nu
        }// mu
        Rcomplex F4gLlgLlgRlgRl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        F4gLlgLlgRlgRl = get_F4(propx1,propy1,propx2,propy2);
        Rcomplex F4pgmuLlgmuLlgnuLlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuLlgmuLlgnuLlgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLlgmuLlgnuLcgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = cpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuLlgmuLlgnuLcgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLcgmuLlgnuLlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = cpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuLcgmuLlgnuLlgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLlgmuLlgRlgLl(0.,0.);
        propy1 = lpropyy;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPL();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F4pgmuLlgmuLlgRlgLl += get_F4(propx1,propx2,propy1,propy2);
        }// mu
        Rcomplex F4pgLlgRlgnuRlgnuRl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glR(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          F4pgLlgRlgnuRlgnuRl += get_F4(propx1,propx2,propy1,propy2);
        }// nu
        Rcomplex F4pgmuLlgmuLlgnuRlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glR(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuLlgmuLlgnuRlgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuRlgmuLlgnuLlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuRlgmuLlgnuLlgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLlgmuLlgnuLlgnuRl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glR(nu);
            F4pgmuLlgmuLlgnuLlgnuRl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLlgmuRlgnuLlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuLlgmuRlgnuLlgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLlgmuLlgRcgLl(0.,0.);
        propy1 = cpropyy;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPL();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F4pgmuLlgmuLlgRcgLl += get_F4(propx1,propx2,propy1,propy2);
        }// mu
        Rcomplex F4pgLcgRlgnuRlgnuRl(0.,0.);
        propx1 = cpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glR(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          F4pgLcgRlgnuRlgnuRl += get_F4(propx1,propx2,propy1,propy2);
        }// nu
        Rcomplex F4pgLlgRlgnuRcgnuRl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = cpropyy;
          propy1 = propy1.glR(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          F4pgLlgRlgnuRcgnuRl += get_F4(propx1,propx2,propy1,propy2);
        }// nu
        Rcomplex F4pgmuLcgmuLlgRlgLl(0.,0.);
        propy1 = lpropyy;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPL();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = cpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F4pgmuLcgmuLlgRlgLl += get_F4(propx1,propx2,propy1,propy2);
        }// mu
        Rcomplex F4pgLlgRlgLlgRl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        propy1 = lpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        F4pgLlgRlgLlgRl = get_F4(propx1,propx2,propy1,propy2);
        Rcomplex F4pgmuRlgmuLlgLlgRl(0.,0.);
        propy1 = lpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F4pgmuRlgmuLlgLlgRl += get_F4(propx1,propx2,propy1,propy2);
        }// mu
        Rcomplex F4pgLlgRlgnuRlgnuLl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glR(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F4pgLlgRlgnuRlgnuLl += get_F4(propx1,propx2,propy1,propy2);
        }// nu
        Rcomplex F4pgmuLlgmuRlgLlgRl(0.,0.);
        propy1 = lpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          F4pgmuLlgmuRlgLlgRl += get_F4(propx1,propx2,propy1,propy2);
        }// mu
        Rcomplex F4pgLlgRlgnuLlgnuRl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          F4pgLlgRlgnuLlgnuRl += get_F4(propx1,propx2,propy1,propy2);
        }// nu
        Rcomplex F4pgLlgRlgLcgRl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        propy1 = cpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        F4pgLlgRlgLcgRl = get_F4(propx1,propx2,propy1,propy2);
        Rcomplex F4pgLcgRlgLlgRl(0.,0.);
        propx1 = cpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        propy1 = lpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        F4pgLcgRlgLlgRl = get_F4(propx1,propx2,propy1,propy2);
        Rcomplex F4pgmuLlgmuRlgnuLcgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = cpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuLlgmuRlgnuLcgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLcgmuLlgnuLlgnuRl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = cpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glR(nu);
            F4pgmuLcgmuLlgnuLlgnuRl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLlgmuRlgnuRlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glR(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuLlgmuRlgnuRlgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuRlgmuLlgnuLlgnuRl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glR(nu);
            F4pgmuRlgmuLlgnuLlgnuRl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLlgmuRlgnuLlgnuRl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glR(nu);
            F4pgmuLlgmuRlgnuLlgnuRl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLlgmuRlgLcgRl(0.,0.);
        propy1 = cpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          F4pgmuLlgmuRlgLcgRl += get_F4(propx1,propx2,propy1,propy2);
        }// mu
        Rcomplex F4pgLcgRlgnuLlgnuRl(0.,0.);
        propx1 = cpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          F4pgLcgRlgnuLlgnuRl += get_F4(propx1,propx2,propy1,propy2);
        }// nu
        Rcomplex F4pgmuRlgmuLlgnuLcgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = cpropyy;
            propy1 = propy1.glL(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuRlgmuLlgnuLcgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuLcgmuLlgnuRlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = cpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glR(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuLcgmuLlgnuRlgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuRlgmuLlgnuRlgnuLl(0.,0.);
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyy;
            propy1 = propy1.glR(nu);
            propy2 = lpropyx;
            propy2 = propy2.glL(nu);
            F4pgmuRlgmuLlgnuRlgnuLl += get_F4(propx1,propx2,propy1,propy2);
          }// nu
        }// mu
        Rcomplex F4pgmuRlgmuLlgLcgRl(0.,0.);
        propy1 = cpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F4pgmuRlgmuLlgLcgRl += get_F4(propx1,propx2,propy1,propy2);
        }// mu
        Rcomplex F4pgLcgRlgnuRlgnuLl(0.,0.);
        propx1 = cpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glR(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F4pgLcgRlgnuRlgnuLl += get_F4(propx1,propx2,propy1,propy2);
        }// nu
        std::vector<Rcomplex> F3xgmuLlgmuLlgnuLl(4, Rcomplex(0.,0.));
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            F3xgmuLlgmuLlgnuLl[nu] += get_F3(propx1,propx2,propy1);
          }// nu
        }// mu
        std::vector<Rcomplex> F3xgmuLcgmuLlgnuLl(4, Rcomplex(0.,0.));
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = cpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            F3xgmuLcgmuLlgnuLl[nu] += get_F3(propx1,propx2,propy1);
          }// nu
        }// mu
        std::vector<Rcomplex> F3xgLlgRlgnuRl(4, Rcomplex(0.,0.));
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyx;
          propy1 = propy1.glR(nu);
          F3xgLlgRlgnuRl[nu] += get_F3(propx1,propx2,propy1);
        }// nu
        std::vector<Rcomplex> F3xgmuRlgmuLlgnuLl(4, Rcomplex(0.,0.));
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            F3xgmuRlgmuLlgnuLl[nu] += get_F3(propx1,propx2,propy1);
          }// nu
        }// mu
        std::vector<Rcomplex> F3xgmuLlgmuRlgnuLl(4, Rcomplex(0.,0.));
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            F3xgmuLlgmuRlgnuLl[nu] += get_F3(propx1,propx2,propy1);
          }// nu
        }// mu
        std::vector<Rcomplex> F3xgLcgRlgnuRl(4, Rcomplex(0.,0.));
        propx1 = cpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyx;
          propy1 = propy1.glR(nu);
          F3xgLcgRlgnuRl[nu] += get_F3(propx1,propx2,propy1);
        }// nu
        Rcomplex F3xgmuLlgmuLlgLl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F3xgmuLlgmuLlgLl += get_F3(propx1,propx2,propy1);
        }// mu
        Rcomplex F3xgmuLlgmuLlgRl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F3xgmuLlgmuLlgRl += get_F3(propx1,propx2,propy1);
        }// mu
        Rcomplex F3xgLlgRlgLl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        F3xgLlgRlgLl = get_F3(propx1,propx2,propy1);
        Rcomplex F3xgLlgRlgRl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        F3xgLlgRlgRl = get_F3(propx1,propx2,propy1);
        Rcomplex F3xgmuLlgmuRlgLl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          F3xgmuLlgmuRlgLl += get_F3(propx1,propx2,propy1);
        }// mu
        Rcomplex F3xgmuLlgmuRlgRl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glR(mu);
          F3xgmuLlgmuRlgRl += get_F3(propx1,propx2,propy1);
        }// mu
        Rcomplex F3xgmuLcgmuLlgRl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = cpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F3xgmuLcgmuLlgRl += get_F3(propx1,propx2,propy1);
        }// mu
        Rcomplex F3xgmuRlgmuLlgRl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F3xgmuRlgmuLlgRl += get_F3(propx1,propx2,propy1);
        }// mu
        Rcomplex F3xgLcgRlgRl(0.,0.);
        propx1 = cpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        F3xgLcgRlgRl = get_F3(propx1,propx2,propy1);
        Rcomplex F3xgmuRlgmuLlgLl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glR(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F3xgmuRlgmuLlgLl += get_F3(propx1,propx2,propy1);
        }// mu
        Rcomplex F3xgmuLcgmuLlgLl(0.,0.);
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = cpropxx;
          propx1 = propx1.glL(mu);
          propx2 = lpropxy;
          propx2 = propx2.glL(mu);
          F3xgmuLcgmuLlgLl += get_F3(propx1,propx2,propy1);
        }// mu
        Rcomplex F3xgLcgRlgLl(0.,0.);
        propx1 = cpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propx2 = lpropxy;
        propx2 = propx2.glPR();
        propx2 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        F3xgLcgRlgLl = get_F3(propx1,propx2,propy1);
        std::vector<Rcomplex> F3ygnuLlgnuLlgmuLl(4, Rcomplex(0.,0.));
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          for ( mu = 0; mu <= 3; ++mu ) {
            propx1 = lpropxy;
            propx1 = propx1.glL(mu);
            F3ygnuLlgnuLlgmuLl[mu] += get_F3(propy1,propy2,propx1);
          }// mu
        }// nu
        std::vector<Rcomplex> F3ygnuLcgnuLlgmuLl(4, Rcomplex(0.,0.));
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = cpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          for ( mu = 0; mu <= 3; ++mu ) {
            propx1 = lpropxy;
            propx1 = propx1.glL(mu);
            F3ygnuLcgnuLlgmuLl[mu] += get_F3(propy1,propy2,propx1);
          }// mu
        }// nu
        std::vector<Rcomplex> F3ygLlgRlgmuRl(4, Rcomplex(0.,0.));
        propy1 = lpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glR(mu);
          F3ygLlgRlgmuRl[mu] += get_F3(propy1,propy2,propx1);
        }// mu
        std::vector<Rcomplex> F3ygnuRlgnuLlgmuLl(4, Rcomplex(0.,0.));
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glR(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          for ( mu = 0; mu <= 3; ++mu ) {
            propx1 = lpropxy;
            propx1 = propx1.glL(mu);
            F3ygnuRlgnuLlgmuLl[mu] += get_F3(propy1,propy2,propx1);
          }// mu
        }// nu
        std::vector<Rcomplex> F3ygnuLlgnuRlgmuLl(4, Rcomplex(0.,0.));
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          for ( mu = 0; mu <= 3; ++mu ) {
            propx1 = lpropxy;
            propx1 = propx1.glL(mu);
            F3ygnuLlgnuRlgmuLl[mu] += get_F3(propy1,propy2,propx1);
          }// mu
        }// nu
        std::vector<Rcomplex> F3ygLcgRlgmuRl(4, Rcomplex(0.,0.));
        propy1 = cpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glR(mu);
          F3ygLcgRlgmuRl[mu] += get_F3(propy1,propy2,propx1);
        }// mu
        Rcomplex F3ygnuLlgnuLlgLl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F3ygnuLlgnuLlgLl += get_F3(propy1,propy2,propx1);
        }// nu
        Rcomplex F3ygnuLlgnuLlgRl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F3ygnuLlgnuLlgRl += get_F3(propy1,propy2,propx1);
        }// nu
        Rcomplex F3ygLlgRlgLl(0.,0.);
        propy1 = lpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        F3ygLlgRlgLl = get_F3(propy1,propy2,propx1);
        Rcomplex F3ygLlgRlgRl(0.,0.);
        propy1 = lpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        F3ygLlgRlgRl = get_F3(propy1,propy2,propx1);
        Rcomplex F3ygnuLlgnuRlgLl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          F3ygnuLlgnuRlgLl += get_F3(propy1,propy2,propx1);
        }// nu
        Rcomplex F3ygnuLlgnuRlgRl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glR(nu);
          F3ygnuLlgnuRlgRl += get_F3(propy1,propy2,propx1);
        }// nu
        Rcomplex F3ygnuLcgnuLlgRl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = cpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F3ygnuLcgnuLlgRl += get_F3(propy1,propy2,propx1);
        }// nu
        Rcomplex F3ygnuRlgnuLlgRl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glR(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F3ygnuRlgnuLlgRl += get_F3(propy1,propy2,propx1);
        }// nu
        Rcomplex F3ygLcgRlgRl(0.,0.);
        propy1 = cpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        F3ygLcgRlgRl = get_F3(propy1,propy2,propx1);
        Rcomplex F3ygnuRlgnuLlgLl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyy;
          propy1 = propy1.glR(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F3ygnuRlgnuLlgLl += get_F3(propy1,propy2,propx1);
        }// nu
        Rcomplex F3ygnuLcgnuLlgLl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = cpropyy;
          propy1 = propy1.glL(nu);
          propy2 = lpropyx;
          propy2 = propy2.glL(nu);
          F3ygnuLcgnuLlgLl += get_F3(propy1,propy2,propx1);
        }// nu
        Rcomplex F3ygLcgRlgLl(0.,0.);
        propy1 = cpropyy;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        propy2 = lpropyx;
        propy2 = propy2.glPR();
        propy2 *= 2.0;
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        F3ygLcgRlgLl = get_F3(propy1,propy2,propx1);
        std::vector<Rcomplex> F2gmuLlgnuLl(16, Rcomplex(0.,0.));
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glL(nu);
            int munu = 4*mu + nu;
            F2gmuLlgnuLl[munu] += Trace(propx1,propy1);
          }// nu
        }// mu
        std::vector<Rcomplex> F2gmuLlgnuRl(16, Rcomplex(0.,0.));
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          for ( nu = 0; nu <= 3; ++nu ) {
            propy1 = lpropyx;
            propy1 = propy1.glR(nu);
            int munu = 4*mu + nu;
            F2gmuLlgnuRl[munu] += Trace(propx1,propy1);
          }// nu
        }// mu
        std::vector<Rcomplex> F2gmuLlgLl(4, Rcomplex(0.,0.));
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          F2gmuLlgLl[mu] += Trace(propx1,propy1);
        }// mu
        std::vector<Rcomplex> F2gmuLlgRl(4, Rcomplex(0.,0.));
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxy;
          propx1 = propx1.glL(mu);
          F2gmuLlgRl[mu] += Trace(propx1,propy1);
        }// mu
        std::vector<Rcomplex> F2gLlgnuLl(4, Rcomplex(0.,0.));
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyx;
          propy1 = propy1.glL(nu);
          F2gLlgnuLl[nu] += Trace(propx1,propy1);
        }// nu
        std::vector<Rcomplex> F2gRlgnuLl(4, Rcomplex(0.,0.));
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        for ( nu = 0; nu <= 3; ++nu ) {
          propy1 = lpropyx;
          propy1 = propy1.glL(nu);
          F2gRlgnuLl[nu] += Trace(propx1,propy1);
        }// nu
        Rcomplex F2gLlgRl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        F2gLlgRl = Trace(propx1,propy1);
        Rcomplex F2gRlgLl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        F2gRlgLl = Trace(propx1,propy1);
        Rcomplex F2gRlgRl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPR();
        propy1 *= 2.0;
        F2gRlgRl = Trace(propx1,propy1);
        Rcomplex F2gLlgLl(0.,0.);
        propx1 = lpropxy;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        propy1 = lpropyx;
        propy1 = propy1.glPL();
        propy1 *= 2.0;
        F2gLlgLl = Trace(propx1,propy1);
        std::vector<Rcomplex> F1xgmuLl(4, Rcomplex(0.,0.));
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = lpropxx;
          propx1 = propx1.glL(mu);
          F1xgmuLl[mu] += Trace(propx1);
        }// mu
        std::vector<Rcomplex> F1xgmuLc(4, Rcomplex(0.,0.));
        for ( mu = 0; mu <= 3; ++mu ) {
          propx1 = cpropxx;
          propx1 = propx1.glL(mu);
          F1xgmuLc[mu] += Trace(propx1);
        }// mu
        Rcomplex F1xgRl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPR();
        propx1 *= 2.0;
        F1xgRl = Trace(propx1);
        Rcomplex F1xgLl(0.,0.);
        propx1 = lpropxx;
        propx1 = propx1.glPL();
        propx1 *= 2.0;
        F1xgLl = Trace(propx1);
	Rcomplex F1xg5l(0.,0.);
	propx1 = lpropxx;
	propx1 = propx1.gl(-5);
	F1xg5l = Trace(propx1);

	for ( mu = 0; mu <= 3; ++mu ) {
	  for ( nu = 0; nu <= 3; ++nu ) {
	    int munu = 4*mu + nu;
	    // __ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_MUNU__
            CON_F2gmuLlgnuLl_F2gmuLlgnuLl[id] += F2gmuLlgnuLl[munu] * F2gmuLlgnuLl[munu];
            CON_F2gmuLlgnuLl_F2gmuLlgnuRl[id] += F2gmuLlgnuLl[munu] * F2gmuLlgnuRl[munu];
            CON_F2gmuLlgnuLl_F2gmuLlgnuRlC[id] += F2gmuLlgnuLl[munu] * conj(F2gmuLlgnuRl[munu]);
            CON_F2gmuLlgnuLl_F2gmuLlgnuLlC[id] += F2gmuLlgnuLl[munu] * conj(F2gmuLlgnuLl[munu]);
            CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLl[id] += F2gmuLlgnuLl[munu] * F1xgmuLl[mu] * F1ygnuLl[nu];
            CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLc[id] += F2gmuLlgnuLl[munu] * F1xgmuLl[mu] * F1ygnuLc[nu];
            CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLl[id] += F2gmuLlgnuLl[munu] * F1xgmuLc[mu] * F1ygnuLl[nu];
            CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLlC[id] += F2gmuLlgnuLl[munu] * F1xgmuLl[mu] * conj(F1ygnuLl[nu]);
            CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLl[id] += F2gmuLlgnuLl[munu] * conj(F1xgmuLl[mu]) * F1ygnuLl[nu];
            CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLcC[id] += F2gmuLlgnuLl[munu] * F1xgmuLl[mu] * conj(F1ygnuLc[nu]);
            CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLl[id] += F2gmuLlgnuLl[munu] * conj(F1xgmuLc[mu]) * F1ygnuLl[nu];
            CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLc[id] += F2gmuLlgnuLl[munu] * conj(F1xgmuLl[mu]) * F1ygnuLc[nu];
            CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLlC[id] += F2gmuLlgnuLl[munu] * F1xgmuLc[mu] * conj(F1ygnuLl[nu]);
            CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLlC[id] += F2gmuLlgnuLl[munu] * conj(F1xgmuLl[mu]) * conj(F1ygnuLl[nu]);
            CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLcC[id] += F2gmuLlgnuLl[munu] * conj(F1xgmuLl[mu]) * conj(F1ygnuLc[nu]);
            CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLlC[id] += F2gmuLlgnuLl[munu] * conj(F1xgmuLc[mu]) * conj(F1ygnuLl[nu]);
	  }
	  // __ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_MU__
          CON_F2gmuLlgLl_F2gmuLlgRl[id] += F2gmuLlgLl[mu] * F2gmuLlgRl[mu];
          CON_F2gLlgnuLl_F2gRlgnuLl[id] += F2gLlgnuLl[mu] * F2gRlgnuLl[mu];
          CON_F2gmuLlgLl_F2gmuLlgRlC[id] += F2gmuLlgLl[mu] * conj(F2gmuLlgRl[mu]);
          CON_F2gLlgnuLl_F2gRlgnuLlC[id] += F2gLlgnuLl[mu] * conj(F2gRlgnuLl[mu]);
          CON_F3xgmuLlgmuLlgnuLl_F1ygnuLl[id] += F3xgmuLlgmuLlgnuLl[mu] * F1ygnuLl[mu];
          CON_F3ygnuLlgnuLlgmuLl_F1xgmuLl[id] += F3ygnuLlgnuLlgmuLl[mu] * F1xgmuLl[mu];
          CON_F3xgmuLcgmuLlgnuLl_F1ygnuLl[id] += F3xgmuLcgmuLlgnuLl[mu] * F1ygnuLl[mu];
          CON_F3ygnuLcgnuLlgmuLl_F1xgmuLl[id] += F3ygnuLcgnuLlgmuLl[mu] * F1xgmuLl[mu];
          CON_F3xgLlgRlgnuRlC_F1ygnuLl[id] += F3xgLlgRlgnuRl[mu] * conj(F1ygnuLl[mu]);
          CON_F3ygLlgRlgmuRlC_F1xgmuLl[id] += F3ygLlgRlgmuRl[mu] * conj(F1xgmuLl[mu]);
          CON_F3xgmuRlgmuLlgnuLl_F1ygnuLl[id] += F3xgmuRlgmuLlgnuLl[mu] * F1ygnuLl[mu];
          CON_F3ygnuRlgnuLlgmuLl_F1xgmuLl[id] += F3ygnuRlgnuLlgmuLl[mu] * F1xgmuLl[mu];
          CON_F3xgmuLlgmuRlgnuLl_F1ygnuLl[id] += F3xgmuLlgmuRlgnuLl[mu] * F1ygnuLl[mu];
          CON_F3ygnuLlgnuRlgmuLl_F1xgmuLl[id] += F3ygnuLlgnuRlgmuLl[mu] * F1xgmuLl[mu];
          CON_F3xgLcgRlgnuRlC_F1ygnuLl[id] += F3xgLcgRlgnuRl[mu] * conj(F1ygnuLl[mu]);
          CON_F3ygLcgRlgmuRlC_F1xgmuLl[id] += F3ygLcgRlgmuRl[mu] * conj(F1xgmuLl[mu]);
          CON_F3xgmuLlgmuLlgnuLl_F1ygnuLc[id] += F3xgmuLlgmuLlgnuLl[mu] * F1ygnuLc[mu];
          CON_F3ygnuLlgnuLlgmuLl_F1xgmuLc[id] += F3ygnuLlgnuLlgmuLl[mu] * F1xgmuLc[mu];
          CON_F3xgmuLlgmuLlgnuLl_F1ygnuLlC[id] += F3xgmuLlgmuLlgnuLl[mu] * conj(F1ygnuLl[mu]);
          CON_F3ygnuLlgnuLlgmuLl_F1xgmuLlC[id] += F3ygnuLlgnuLlgmuLl[mu] * conj(F1xgmuLl[mu]);
          CON_F3xgmuLlgmuLlgnuLl_F1ygnuLcC[id] += F3xgmuLlgmuLlgnuLl[mu] * conj(F1ygnuLc[mu]);
          CON_F3ygnuLlgnuLlgmuLl_F1xgmuLcC[id] += F3ygnuLlgnuLlgmuLl[mu] * conj(F1xgmuLc[mu]);
          CON_F3xgmuLcgmuLlgnuLl_F1ygnuLlC[id] += F3xgmuLcgmuLlgnuLl[mu] * conj(F1ygnuLl[mu]);
          CON_F3ygnuLcgnuLlgmuLl_F1xgmuLlC[id] += F3ygnuLcgnuLlgmuLl[mu] * conj(F1xgmuLl[mu]);
          CON_F3xgLlgRlgnuRlC_F1ygnuLlC[id] += F3xgLlgRlgnuRl[mu] * F1ygnuLl[mu];
          CON_F3ygLlgRlgmuRlC_F1xgmuLlC[id] += F3ygLlgRlgmuRl[mu] * F1xgmuLl[mu];
          CON_F3xgmuRlgmuLlgnuLl_F1ygnuLlC[id] += F3xgmuRlgmuLlgnuLl[mu] * conj(F1ygnuLl[mu]);
          CON_F3ygnuRlgnuLlgmuLl_F1xgmuLlC[id] += F3ygnuRlgnuLlgmuLl[mu] * conj(F1xgmuLl[mu]);
          CON_F3xgmuLlgmuRlgnuLl_F1ygnuLlC[id] += F3xgmuLlgmuRlgnuLl[mu] * conj(F1ygnuLl[mu]);
          CON_F3ygnuLlgnuRlgmuLl_F1xgmuLlC[id] += F3ygnuLlgnuRlgmuLl[mu] * conj(F1xgmuLl[mu]);
          CON_F3xgLcgRlgnuRlC_F1ygnuLlC[id] += F3xgLcgRlgnuRl[mu] * F1ygnuLl[mu];
          CON_F3ygLcgRlgmuRlC_F1xgmuLlC[id] += F3ygLcgRlgmuRl[mu] * F1xgmuLl[mu];
          CON_F3xgLlgRlgnuRlC_F1ygnuLc[id] += F3xgLlgRlgnuRl[mu] * conj(F1ygnuLc[mu]);
          CON_F3ygLlgRlgmuRlC_F1xgmuLc[id] += F3ygLlgRlgmuRl[mu] * conj(F1xgmuLc[mu]);
          CON_F3xgLlgRlgnuRlC_F1ygnuLcC[id] += F3xgLlgRlgnuRl[mu] * F1ygnuLc[mu];
          CON_F3ygLlgRlgmuRlC_F1xgmuLcC[id] += F3ygLlgRlgmuRl[mu] * F1xgmuLc[mu];
          CON_F3xgmuLlgmuRlgnuLl_F1ygnuLc[id] += F3xgmuLlgmuRlgnuLl[mu] * F1ygnuLc[mu];
          CON_F3ygnuLlgnuRlgmuLl_F1xgmuLc[id] += F3ygnuLlgnuRlgmuLl[mu] * F1xgmuLc[mu];
          CON_F3xgmuLlgmuRlgnuLl_F1ygnuLcC[id] += F3xgmuLlgmuRlgnuLl[mu] * conj(F1ygnuLc[mu]);
          CON_F3ygnuLlgnuRlgmuLl_F1xgmuLcC[id] += F3ygnuLlgnuRlgmuLl[mu] * conj(F1xgmuLc[mu]);
          CON_F3xgmuRlgmuLlgnuLl_F1ygnuLc[id] += F3xgmuRlgmuLlgnuLl[mu] * F1ygnuLc[mu];
          CON_F3ygnuRlgnuLlgmuLl_F1xgmuLc[id] += F3ygnuRlgnuLlgmuLl[mu] * F1xgmuLc[mu];
          CON_F3xgmuRlgmuLlgnuLl_F1ygnuLcC[id] += F3xgmuRlgmuLlgnuLl[mu] * conj(F1ygnuLc[mu]);
          CON_F3ygnuRlgnuLlgmuLl_F1xgmuLcC[id] += F3ygnuRlgnuLlgmuLl[mu] * conj(F1xgmuLc[mu]);
          CON_F2gmuLlgLl_F1xgmuLl_F1ygRl[id] += F2gmuLlgLl[mu] * F1xgmuLl[mu] * F1ygRl;
          CON_F2gLlgnuLl_F1xgRl_F1ygnuLl[id] += F2gLlgnuLl[mu] * F1xgRl * F1ygnuLl[mu];
          CON_F2gmuLlgRl_F1xgmuLl_F1ygLl[id] += F2gmuLlgRl[mu] * F1xgmuLl[mu] * F1ygLl;
          CON_F2gRlgnuLl_F1xgLl_F1ygnuLl[id] += F2gRlgnuLl[mu] * F1xgLl * F1ygnuLl[mu];
          CON_F2gmuLlgLl_F1xgmuLlC_F1ygRl[id] += F2gmuLlgLl[mu] * conj(F1xgmuLl[mu]) * F1ygRl;
          CON_F2gLlgnuLl_F1xgRl_F1ygnuLlC[id] += F2gLlgnuLl[mu] * F1xgRl * conj(F1ygnuLl[mu]);
          CON_F2gmuLlgRl_F1xgmuLlC_F1ygLl[id] += F2gmuLlgRl[mu] * conj(F1xgmuLl[mu]) * F1ygLl;
          CON_F2gRlgnuLl_F1xgLl_F1ygnuLlC[id] += F2gRlgnuLl[mu] * F1xgLl * conj(F1ygnuLl[mu]);
          CON_F2gmuLlgRl_F1xgmuLc_F1ygLl[id] += F2gmuLlgRl[mu] * F1xgmuLc[mu] * F1ygLl;
          CON_F2gRlgnuLl_F1xgLl_F1ygnuLc[id] += F2gRlgnuLl[mu] * F1xgLl * F1ygnuLc[mu];
          CON_F2gmuLlgRl_F1xgmuLcC_F1ygLl[id] += F2gmuLlgRl[mu] * conj(F1xgmuLc[mu]) * F1ygLl;
          CON_F2gRlgnuLl_F1xgLl_F1ygnuLcC[id] += F2gRlgnuLl[mu] * F1xgLl * conj(F1ygnuLc[mu]);
          CON_F2gmuLlgLl_F1xgmuLc_F1ygRl[id] += F2gmuLlgLl[mu] * F1xgmuLc[mu] * F1ygRl;
          CON_F2gLlgnuLl_F1xgRl_F1ygnuLc[id] += F2gLlgnuLl[mu] * F1xgRl * F1ygnuLc[mu];
          CON_F2gmuLlgLl_F1xgmuLcC_F1ygRl[id] += F2gmuLlgLl[mu] * conj(F1xgmuLc[mu]) * F1ygRl;
          CON_F2gLlgnuLl_F1xgRl_F1ygnuLcC[id] += F2gLlgnuLl[mu] * F1xgRl * conj(F1ygnuLc[mu]);
          CON_F2gmuLlgRl_F1xgmuLl[id] += F2gmuLlgRl[mu] * F1xgmuLl[mu];
          CON_F2gRlgnuLl_F1ygnuLl[id] += F2gRlgnuLl[mu] * F1ygnuLl[mu];
          CON_F2gmuLlgLl_F1xgmuLl[id] += F2gmuLlgLl[mu] * F1xgmuLl[mu];
          CON_F2gLlgnuLl_F1ygnuLl[id] += F2gLlgnuLl[mu] * F1ygnuLl[mu];
          CON_F2gmuLlgRl_F1xgmuLlC[id] += F2gmuLlgRl[mu] * conj(F1xgmuLl[mu]);
          CON_F2gRlgnuLl_F1ygnuLlC[id] += F2gRlgnuLl[mu] * conj(F1ygnuLl[mu]);
          CON_F2gmuLlgLl_F1xgmuLlC[id] += F2gmuLlgLl[mu] * conj(F1xgmuLl[mu]);
          CON_F2gLlgnuLl_F1ygnuLlC[id] += F2gLlgnuLl[mu] * conj(F1ygnuLl[mu]);
          CON_F2gmuLlgLl_F1xgmuLc[id] += F2gmuLlgLl[mu] * F1xgmuLc[mu];
          CON_F2gLlgnuLl_F1ygnuLc[id] += F2gLlgnuLl[mu] * F1ygnuLc[mu];
          CON_F2gmuLlgLl_F1xgmuLcC[id] += F2gmuLlgLl[mu] * conj(F1xgmuLc[mu]);
          CON_F2gLlgnuLl_F1ygnuLcC[id] += F2gLlgnuLl[mu] * conj(F1ygnuLc[mu]);
          CON_F2gmuLlgRl_F1xgmuLc[id] += F2gmuLlgRl[mu] * F1xgmuLc[mu];
          CON_F2gRlgnuLl_F1ygnuLc[id] += F2gRlgnuLl[mu] * F1ygnuLc[mu];
          CON_F2gmuLlgRl_F1xgmuLcC[id] += F2gmuLlgRl[mu] * conj(F1xgmuLc[mu]);
          CON_F2gRlgnuLl_F1ygnuLcC[id] += F2gRlgnuLl[mu] * conj(F1ygnuLc[mu]);
	}
	// __ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_SCL__
        CON_F4gmuLlgnuLlgmuLlgnuLl[id] += F4gmuLlgnuLlgmuLlgnuLl;
        CON_F4gmuLlgLlgmuLlgRl[id] += F4gmuLlgLlgmuLlgRl;
        CON_F4gLlgnuLlgRlgnuLl[id] += F4gLlgnuLlgRlgnuLl;
        CON_F4gmuLlgnuLlgmuLlgnuRl[id] += F4gmuLlgnuLlgmuLlgnuRl;
        CON_F4gmuLlgnuLlgmuRlgnuLl[id] += F4gmuLlgnuLlgmuRlgnuLl;
        CON_F4gmuLlgRlgmuRlgLl[id] += F4gmuLlgRlgmuRlgLl;
        CON_F4gLlgnuLlgRlgnuRl[id] += F4gLlgnuLlgRlgnuRl;
        CON_F4gmuLlgnuLlgmuRlgnuRl[id] += F4gmuLlgnuLlgmuRlgnuRl;
        CON_F4gLlgLlgRlgRl[id] += F4gLlgLlgRlgRl;
        CON_F2gLlgRl_F2gRlgLl[id] += F2gLlgRl * F2gRlgLl;
        CON_F4pgmuLlgmuLlgnuLlgnuLl[id] += F4pgmuLlgmuLlgnuLlgnuLl;
        CON_F4pgmuLlgmuLlgnuLcgnuLl[id] += F4pgmuLlgmuLlgnuLcgnuLl;
        CON_F4pgmuLcgmuLlgnuLlgnuLl[id] += F4pgmuLcgmuLlgnuLlgnuLl;
        CON_F4pgmuLlgmuLlgRlgLl[id] += F4pgmuLlgmuLlgRlgLl;
        CON_F4pgLlgRlgnuRlgnuRl[id] += F4pgLlgRlgnuRlgnuRl;
        CON_F4pgmuLlgmuLlgnuRlgnuLl[id] += F4pgmuLlgmuLlgnuRlgnuLl;
        CON_F4pgmuRlgmuLlgnuLlgnuLl[id] += F4pgmuRlgmuLlgnuLlgnuLl;
        CON_F4pgmuLlgmuLlgnuLlgnuRl[id] += F4pgmuLlgmuLlgnuLlgnuRl;
        CON_F4pgmuLlgmuRlgnuLlgnuLl[id] += F4pgmuLlgmuRlgnuLlgnuLl;
        CON_F4pgmuLlgmuLlgRcgLl[id] += F4pgmuLlgmuLlgRcgLl;
        CON_F4pgLcgRlgnuRlgnuRl[id] += F4pgLcgRlgnuRlgnuRl;
        CON_F4pgLlgRlgnuRcgnuRl[id] += F4pgLlgRlgnuRcgnuRl;
        CON_F4pgmuLcgmuLlgRlgLl[id] += F4pgmuLcgmuLlgRlgLl;
        CON_F4pgLlgRlgLlgRl[id] += F4pgLlgRlgLlgRl;
        CON_F4pgmuRlgmuLlgLlgRl[id] += F4pgmuRlgmuLlgLlgRl;
        CON_F4pgLlgRlgnuRlgnuLl[id] += F4pgLlgRlgnuRlgnuLl;
        CON_F4pgmuLlgmuRlgLlgRl[id] += F4pgmuLlgmuRlgLlgRl;
        CON_F4pgLlgRlgnuLlgnuRl[id] += F4pgLlgRlgnuLlgnuRl;
        CON_F4pgLlgRlgLcgRl[id] += F4pgLlgRlgLcgRl;
        CON_F4pgLcgRlgLlgRl[id] += F4pgLcgRlgLlgRl;
        CON_F4pgmuLlgmuRlgnuLcgnuLl[id] += F4pgmuLlgmuRlgnuLcgnuLl;
        CON_F4pgmuLcgmuLlgnuLlgnuRl[id] += F4pgmuLcgmuLlgnuLlgnuRl;
        CON_F4pgmuLlgmuRlgnuRlgnuLl[id] += F4pgmuLlgmuRlgnuRlgnuLl;
        CON_F4pgmuRlgmuLlgnuLlgnuRl[id] += F4pgmuRlgmuLlgnuLlgnuRl;
        CON_F4pgmuLlgmuRlgnuLlgnuRl[id] += F4pgmuLlgmuRlgnuLlgnuRl;
        CON_F4pgmuLlgmuRlgLcgRl[id] += F4pgmuLlgmuRlgLcgRl;
        CON_F4pgLcgRlgnuLlgnuRl[id] += F4pgLcgRlgnuLlgnuRl;
        CON_F4pgmuRlgmuLlgnuLcgnuLl[id] += F4pgmuRlgmuLlgnuLcgnuLl;
        CON_F4pgmuLcgmuLlgnuRlgnuLl[id] += F4pgmuLcgmuLlgnuRlgnuLl;
        CON_F4pgmuRlgmuLlgnuRlgnuLl[id] += F4pgmuRlgmuLlgnuRlgnuLl;
        CON_F4pgmuRlgmuLlgLcgRl[id] += F4pgmuRlgmuLlgLcgRl;
        CON_F4pgLcgRlgnuRlgnuLl[id] += F4pgLcgRlgnuRlgnuLl;
        CON_F3xgmuLlgmuLlgLl_F1ygRl[id] += F3xgmuLlgmuLlgLl * F1ygRl;
        CON_F3ygnuLlgnuLlgLl_F1xgRl[id] += F3ygnuLlgnuLlgLl * F1xgRl;
        CON_F3xgmuLlgmuLlgRl_F1ygLl[id] += F3xgmuLlgmuLlgRl * F1ygLl;
        CON_F3ygnuLlgnuLlgRl_F1xgLl[id] += F3ygnuLlgnuLlgRl * F1xgLl;
        CON_F3xgLlgRlgLlC_F1ygRl[id] += F3xgLlgRlgLl * conj(F1ygRl);
        CON_F3ygLlgRlgLlC_F1xgRl[id] += F3ygLlgRlgLl * conj(F1xgRl);
        CON_F3xgLlgRlgRlC_F1ygLl[id] += F3xgLlgRlgRl * conj(F1ygLl);
        CON_F3ygLlgRlgRlC_F1xgLl[id] += F3ygLlgRlgRl * conj(F1xgLl);
        CON_F3xgmuLlgmuRlgLl_F1ygRl[id] += F3xgmuLlgmuRlgLl * F1ygRl;
        CON_F3ygnuLlgnuRlgLl_F1xgRl[id] += F3ygnuLlgnuRlgLl * F1xgRl;
        CON_F3xgmuLlgmuRlgRl_F1ygLl[id] += F3xgmuLlgmuRlgRl * F1ygLl;
        CON_F3ygnuLlgnuRlgRl_F1xgLl[id] += F3ygnuLlgnuRlgRl * F1xgLl;
        CON_F3xgmuLcgmuLlgRl_F1ygLl[id] += F3xgmuLcgmuLlgRl * F1ygLl;
        CON_F3ygnuLcgnuLlgRl_F1xgLl[id] += F3ygnuLcgnuLlgRl * F1xgLl;
        CON_F3xgmuRlgmuLlgRl_F1ygLl[id] += F3xgmuRlgmuLlgRl * F1ygLl;
        CON_F3ygnuRlgnuLlgRl_F1xgLl[id] += F3ygnuRlgnuLlgRl * F1xgLl;
        CON_F3xgLcgRlgRlC_F1ygLl[id] += F3xgLcgRlgRl * conj(F1ygLl);
        CON_F3ygLcgRlgRlC_F1xgLl[id] += F3ygLcgRlgRl * conj(F1xgLl);
        CON_F3xgmuRlgmuLlgLl_F1ygRl[id] += F3xgmuRlgmuLlgLl * F1ygRl;
        CON_F3ygnuRlgnuLlgLl_F1xgRl[id] += F3ygnuRlgnuLlgLl * F1xgRl;
        CON_F3xgmuLcgmuLlgLl_F1ygRl[id] += F3xgmuLcgmuLlgLl * F1ygRl;
        CON_F3ygnuLcgnuLlgLl_F1xgRl[id] += F3ygnuLcgnuLlgLl * F1xgRl;
        CON_F3xgLcgRlgLlC_F1ygRl[id] += F3xgLcgRlgLl * conj(F1ygRl);
        CON_F3ygLcgRlgLlC_F1xgRl[id] += F3ygLcgRlgLl * conj(F1xgRl);
        CON_F2gLlgRl_F1xgRl_F1ygLl[id] += F2gLlgRl * F1xgRl * F1ygLl;
        CON_F2gRlgLl_F1xgLl_F1ygRl[id] += F2gRlgLl * F1xgLl * F1ygRl;
        CON_F2gRlgRl_F1xgLl_F1ygLl[id] += F2gRlgRl * F1xgLl * F1ygLl;
        CON_F2gLlgLl_F1xgRl_F1ygRl[id] += F2gLlgLl * F1xgRl * F1ygRl;
        CON_F3xgmuLlgmuLlgRl[id] += F3xgmuLlgmuLlgRl;
        CON_F3ygnuLlgnuLlgRl[id] += F3ygnuLlgnuLlgRl;
        CON_F3xgmuLlgmuLlgLl[id] += F3xgmuLlgmuLlgLl;
        CON_F3ygnuLlgnuLlgLl[id] += F3ygnuLlgnuLlgLl;
        CON_F3xgLlgRlgRl[id] += F3xgLlgRlgRl;
        CON_F3ygLlgRlgRl[id] += F3ygLlgRlgRl;
        CON_F3xgLlgRlgLl[id] += F3xgLlgRlgLl;
        CON_F3ygLlgRlgLl[id] += F3ygLlgRlgLl;
        CON_F3xgmuLlgmuRlgRl[id] += F3xgmuLlgmuRlgRl;
        CON_F3ygnuLlgnuRlgRl[id] += F3ygnuLlgnuRlgRl;
        CON_F3xgmuLlgmuRlgLl[id] += F3xgmuLlgmuRlgLl;
        CON_F3ygnuLlgnuRlgLl[id] += F3ygnuLlgnuRlgLl;
        CON_F3xgmuRlgmuLlgRl[id] += F3xgmuRlgmuLlgRl;
        CON_F3ygnuRlgnuLlgRl[id] += F3ygnuRlgnuLlgRl;
        CON_F3xgmuRlgmuLlgLl[id] += F3xgmuRlgmuLlgLl;
        CON_F3ygnuRlgnuLlgLl[id] += F3ygnuRlgnuLlgLl;
        CON_F3xgmuLcgmuLlgLl[id] += F3xgmuLcgmuLlgLl;
        CON_F3ygnuLcgnuLlgLl[id] += F3ygnuLcgnuLlgLl;
        CON_F3xgLcgRlgLl[id] += F3xgLcgRlgLl;
        CON_F3ygLcgRlgLl[id] += F3ygLcgRlgLl;
        CON_F3xgmuLcgmuLlgRl[id] += F3xgmuLcgmuLlgRl;
        CON_F3ygnuLcgnuLlgRl[id] += F3ygnuLcgnuLlgRl;
        CON_F3xgLcgRlgRl[id] += F3xgLcgRlgRl;
        CON_F3ygLcgRlgRl[id] += F3ygLcgRlgRl;
        CON_F2gRlgRl_F1xgLl[id] += F2gRlgRl * F1xgLl;
        CON_F2gRlgRl_F1ygLl[id] += F2gRlgRl * F1ygLl;
        CON_F2gRlgLl_F1xgLl[id] += F2gRlgLl * F1xgLl;
        CON_F2gLlgRl_F1ygLl[id] += F2gLlgRl * F1ygLl;
        CON_F2gLlgRl_F1xgRl[id] += F2gLlgRl * F1xgRl;
        CON_F2gRlgLl_F1ygRl[id] += F2gRlgLl * F1ygRl;
        CON_F2gLlgLl_F1xgRl[id] += F2gLlgLl * F1xgRl;
        CON_F2gLlgLl_F1ygRl[id] += F2gLlgLl * F1ygRl;
        CON_F2gLlgRl[id] += F2gLlgRl;
        CON_F2gRlgLl[id] += F2gRlgLl;
        CON_F2gLlgLl[id] += F2gLlgLl;
        CON_F2gRlgRl[id] += F2gRlgRl;
	corr_disc[rel_t] += F1xg5l * F1yg5l;
	s.nextSite();
      }// loop over sites
    }// loop over source points
    // __SUM_OVER_ALL_NODES__
    glb_sum(CON_F4gmuLlgnuLlgmuLlgnuLl.data(), CON_F4gmuLlgnuLlgmuLlgnuLl.size());
    glb_sum(CON_F4gmuLlgLlgmuLlgRl.data(), CON_F4gmuLlgLlgmuLlgRl.size());
    glb_sum(CON_F4gLlgnuLlgRlgnuLl.data(), CON_F4gLlgnuLlgRlgnuLl.size());
    glb_sum(CON_F4gmuLlgnuLlgmuLlgnuRl.data(), CON_F4gmuLlgnuLlgmuLlgnuRl.size());
    glb_sum(CON_F4gmuLlgnuLlgmuRlgnuLl.data(), CON_F4gmuLlgnuLlgmuRlgnuLl.size());
    glb_sum(CON_F4gmuLlgRlgmuRlgLl.data(), CON_F4gmuLlgRlgmuRlgLl.size());
    glb_sum(CON_F4gLlgnuLlgRlgnuRl.data(), CON_F4gLlgnuLlgRlgnuRl.size());
    glb_sum(CON_F4gmuLlgnuLlgmuRlgnuRl.data(), CON_F4gmuLlgnuLlgmuRlgnuRl.size());
    glb_sum(CON_F4gLlgLlgRlgRl.data(), CON_F4gLlgLlgRlgRl.size());
    glb_sum(CON_F2gmuLlgnuLl_F2gmuLlgnuLl.data(), CON_F2gmuLlgnuLl_F2gmuLlgnuLl.size());
    glb_sum(CON_F2gmuLlgnuLl_F2gmuLlgnuRl.data(), CON_F2gmuLlgnuLl_F2gmuLlgnuRl.size());
    glb_sum(CON_F2gmuLlgnuLl_F2gmuLlgnuRlC.data(), CON_F2gmuLlgnuLl_F2gmuLlgnuRlC.size());
    glb_sum(CON_F2gmuLlgLl_F2gmuLlgRl.data(), CON_F2gmuLlgLl_F2gmuLlgRl.size());
    glb_sum(CON_F2gLlgnuLl_F2gRlgnuLl.data(), CON_F2gLlgnuLl_F2gRlgnuLl.size());
    glb_sum(CON_F2gmuLlgnuLl_F2gmuLlgnuLlC.data(), CON_F2gmuLlgnuLl_F2gmuLlgnuLlC.size());
    glb_sum(CON_F2gLlgRl_F2gRlgLl.data(), CON_F2gLlgRl_F2gRlgLl.size());
    glb_sum(CON_F2gmuLlgLl_F2gmuLlgRlC.data(), CON_F2gmuLlgLl_F2gmuLlgRlC.size());
    glb_sum(CON_F2gLlgnuLl_F2gRlgnuLlC.data(), CON_F2gLlgnuLl_F2gRlgnuLlC.size());
    glb_sum(CON_F4pgmuLlgmuLlgnuLlgnuLl.data(), CON_F4pgmuLlgmuLlgnuLlgnuLl.size());
    glb_sum(CON_F4pgmuLlgmuLlgnuLcgnuLl.data(), CON_F4pgmuLlgmuLlgnuLcgnuLl.size());
    glb_sum(CON_F4pgmuLcgmuLlgnuLlgnuLl.data(), CON_F4pgmuLcgmuLlgnuLlgnuLl.size());
    glb_sum(CON_F4pgmuLlgmuLlgRlgLl.data(), CON_F4pgmuLlgmuLlgRlgLl.size());
    glb_sum(CON_F4pgLlgRlgnuRlgnuRl.data(), CON_F4pgLlgRlgnuRlgnuRl.size());
    glb_sum(CON_F4pgmuLlgmuLlgnuRlgnuLl.data(), CON_F4pgmuLlgmuLlgnuRlgnuLl.size());
    glb_sum(CON_F4pgmuRlgmuLlgnuLlgnuLl.data(), CON_F4pgmuRlgmuLlgnuLlgnuLl.size());
    glb_sum(CON_F4pgmuLlgmuLlgnuLlgnuRl.data(), CON_F4pgmuLlgmuLlgnuLlgnuRl.size());
    glb_sum(CON_F4pgmuLlgmuRlgnuLlgnuLl.data(), CON_F4pgmuLlgmuRlgnuLlgnuLl.size());
    glb_sum(CON_F4pgmuLlgmuLlgRcgLl.data(), CON_F4pgmuLlgmuLlgRcgLl.size());
    glb_sum(CON_F4pgLcgRlgnuRlgnuRl.data(), CON_F4pgLcgRlgnuRlgnuRl.size());
    glb_sum(CON_F4pgLlgRlgnuRcgnuRl.data(), CON_F4pgLlgRlgnuRcgnuRl.size());
    glb_sum(CON_F4pgmuLcgmuLlgRlgLl.data(), CON_F4pgmuLcgmuLlgRlgLl.size());
    glb_sum(CON_F4pgLlgRlgLlgRl.data(), CON_F4pgLlgRlgLlgRl.size());
    glb_sum(CON_F4pgmuRlgmuLlgLlgRl.data(), CON_F4pgmuRlgmuLlgLlgRl.size());
    glb_sum(CON_F4pgLlgRlgnuRlgnuLl.data(), CON_F4pgLlgRlgnuRlgnuLl.size());
    glb_sum(CON_F4pgmuLlgmuRlgLlgRl.data(), CON_F4pgmuLlgmuRlgLlgRl.size());
    glb_sum(CON_F4pgLlgRlgnuLlgnuRl.data(), CON_F4pgLlgRlgnuLlgnuRl.size());
    glb_sum(CON_F4pgLlgRlgLcgRl.data(), CON_F4pgLlgRlgLcgRl.size());
    glb_sum(CON_F4pgLcgRlgLlgRl.data(), CON_F4pgLcgRlgLlgRl.size());
    glb_sum(CON_F4pgmuLlgmuRlgnuLcgnuLl.data(), CON_F4pgmuLlgmuRlgnuLcgnuLl.size());
    glb_sum(CON_F4pgmuLcgmuLlgnuLlgnuRl.data(), CON_F4pgmuLcgmuLlgnuLlgnuRl.size());
    glb_sum(CON_F4pgmuLlgmuRlgnuRlgnuLl.data(), CON_F4pgmuLlgmuRlgnuRlgnuLl.size());
    glb_sum(CON_F4pgmuRlgmuLlgnuLlgnuRl.data(), CON_F4pgmuRlgmuLlgnuLlgnuRl.size());
    glb_sum(CON_F4pgmuLlgmuRlgnuLlgnuRl.data(), CON_F4pgmuLlgmuRlgnuLlgnuRl.size());
    glb_sum(CON_F4pgmuLlgmuRlgLcgRl.data(), CON_F4pgmuLlgmuRlgLcgRl.size());
    glb_sum(CON_F4pgLcgRlgnuLlgnuRl.data(), CON_F4pgLcgRlgnuLlgnuRl.size());
    glb_sum(CON_F4pgmuRlgmuLlgnuLcgnuLl.data(), CON_F4pgmuRlgmuLlgnuLcgnuLl.size());
    glb_sum(CON_F4pgmuLcgmuLlgnuRlgnuLl.data(), CON_F4pgmuLcgmuLlgnuRlgnuLl.size());
    glb_sum(CON_F4pgmuRlgmuLlgnuRlgnuLl.data(), CON_F4pgmuRlgmuLlgnuRlgnuLl.size());
    glb_sum(CON_F4pgmuRlgmuLlgLcgRl.data(), CON_F4pgmuRlgmuLlgLcgRl.size());
    glb_sum(CON_F4pgLcgRlgnuRlgnuLl.data(), CON_F4pgLcgRlgnuRlgnuLl.size());
    glb_sum(CON_F3xgmuLlgmuLlgnuLl_F1ygnuLl.data(), CON_F3xgmuLlgmuLlgnuLl_F1ygnuLl.size());
    glb_sum(CON_F3ygnuLlgnuLlgmuLl_F1xgmuLl.data(), CON_F3ygnuLlgnuLlgmuLl_F1xgmuLl.size());
    glb_sum(CON_F3xgmuLcgmuLlgnuLl_F1ygnuLl.data(), CON_F3xgmuLcgmuLlgnuLl_F1ygnuLl.size());
    glb_sum(CON_F3ygnuLcgnuLlgmuLl_F1xgmuLl.data(), CON_F3ygnuLcgnuLlgmuLl_F1xgmuLl.size());
    glb_sum(CON_F3xgLlgRlgnuRlC_F1ygnuLl.data(), CON_F3xgLlgRlgnuRlC_F1ygnuLl.size());
    glb_sum(CON_F3ygLlgRlgmuRlC_F1xgmuLl.data(), CON_F3ygLlgRlgmuRlC_F1xgmuLl.size());
    glb_sum(CON_F3xgmuRlgmuLlgnuLl_F1ygnuLl.data(), CON_F3xgmuRlgmuLlgnuLl_F1ygnuLl.size());
    glb_sum(CON_F3ygnuRlgnuLlgmuLl_F1xgmuLl.data(), CON_F3ygnuRlgnuLlgmuLl_F1xgmuLl.size());
    glb_sum(CON_F3xgmuLlgmuRlgnuLl_F1ygnuLl.data(), CON_F3xgmuLlgmuRlgnuLl_F1ygnuLl.size());
    glb_sum(CON_F3ygnuLlgnuRlgmuLl_F1xgmuLl.data(), CON_F3ygnuLlgnuRlgmuLl_F1xgmuLl.size());
    glb_sum(CON_F3xgLcgRlgnuRlC_F1ygnuLl.data(), CON_F3xgLcgRlgnuRlC_F1ygnuLl.size());
    glb_sum(CON_F3ygLcgRlgmuRlC_F1xgmuLl.data(), CON_F3ygLcgRlgmuRlC_F1xgmuLl.size());
    glb_sum(CON_F3xgmuLlgmuLlgnuLl_F1ygnuLc.data(), CON_F3xgmuLlgmuLlgnuLl_F1ygnuLc.size());
    glb_sum(CON_F3ygnuLlgnuLlgmuLl_F1xgmuLc.data(), CON_F3ygnuLlgnuLlgmuLl_F1xgmuLc.size());
    glb_sum(CON_F3xgmuLlgmuLlgnuLl_F1ygnuLlC.data(), CON_F3xgmuLlgmuLlgnuLl_F1ygnuLlC.size());
    glb_sum(CON_F3ygnuLlgnuLlgmuLl_F1xgmuLlC.data(), CON_F3ygnuLlgnuLlgmuLl_F1xgmuLlC.size());
    glb_sum(CON_F3xgmuLlgmuLlgLl_F1ygRl.data(), CON_F3xgmuLlgmuLlgLl_F1ygRl.size());
    glb_sum(CON_F3ygnuLlgnuLlgLl_F1xgRl.data(), CON_F3ygnuLlgnuLlgLl_F1xgRl.size());
    glb_sum(CON_F3xgmuLlgmuLlgRl_F1ygLl.data(), CON_F3xgmuLlgmuLlgRl_F1ygLl.size());
    glb_sum(CON_F3ygnuLlgnuLlgRl_F1xgLl.data(), CON_F3ygnuLlgnuLlgRl_F1xgLl.size());
    glb_sum(CON_F3xgmuLlgmuLlgnuLl_F1ygnuLcC.data(), CON_F3xgmuLlgmuLlgnuLl_F1ygnuLcC.size());
    glb_sum(CON_F3ygnuLlgnuLlgmuLl_F1xgmuLcC.data(), CON_F3ygnuLlgnuLlgmuLl_F1xgmuLcC.size());
    glb_sum(CON_F3xgmuLcgmuLlgnuLl_F1ygnuLlC.data(), CON_F3xgmuLcgmuLlgnuLl_F1ygnuLlC.size());
    glb_sum(CON_F3ygnuLcgnuLlgmuLl_F1xgmuLlC.data(), CON_F3ygnuLcgnuLlgmuLl_F1xgmuLlC.size());
    glb_sum(CON_F3xgLlgRlgnuRlC_F1ygnuLlC.data(), CON_F3xgLlgRlgnuRlC_F1ygnuLlC.size());
    glb_sum(CON_F3ygLlgRlgmuRlC_F1xgmuLlC.data(), CON_F3ygLlgRlgmuRlC_F1xgmuLlC.size());
    glb_sum(CON_F3xgmuRlgmuLlgnuLl_F1ygnuLlC.data(), CON_F3xgmuRlgmuLlgnuLl_F1ygnuLlC.size());
    glb_sum(CON_F3ygnuRlgnuLlgmuLl_F1xgmuLlC.data(), CON_F3ygnuRlgnuLlgmuLl_F1xgmuLlC.size());
    glb_sum(CON_F3xgmuLlgmuRlgnuLl_F1ygnuLlC.data(), CON_F3xgmuLlgmuRlgnuLl_F1ygnuLlC.size());
    glb_sum(CON_F3ygnuLlgnuRlgmuLl_F1xgmuLlC.data(), CON_F3ygnuLlgnuRlgmuLl_F1xgmuLlC.size());
    glb_sum(CON_F3xgLcgRlgnuRlC_F1ygnuLlC.data(), CON_F3xgLcgRlgnuRlC_F1ygnuLlC.size());
    glb_sum(CON_F3ygLcgRlgmuRlC_F1xgmuLlC.data(), CON_F3ygLcgRlgmuRlC_F1xgmuLlC.size());
    glb_sum(CON_F3xgLlgRlgnuRlC_F1ygnuLc.data(), CON_F3xgLlgRlgnuRlC_F1ygnuLc.size());
    glb_sum(CON_F3ygLlgRlgmuRlC_F1xgmuLc.data(), CON_F3ygLlgRlgmuRlC_F1xgmuLc.size());
    glb_sum(CON_F3xgLlgRlgLlC_F1ygRl.data(), CON_F3xgLlgRlgLlC_F1ygRl.size());
    glb_sum(CON_F3ygLlgRlgLlC_F1xgRl.data(), CON_F3ygLlgRlgLlC_F1xgRl.size());
    glb_sum(CON_F3xgLlgRlgRlC_F1ygLl.data(), CON_F3xgLlgRlgRlC_F1ygLl.size());
    glb_sum(CON_F3ygLlgRlgRlC_F1xgLl.data(), CON_F3ygLlgRlgRlC_F1xgLl.size());
    glb_sum(CON_F3xgLlgRlgnuRlC_F1ygnuLcC.data(), CON_F3xgLlgRlgnuRlC_F1ygnuLcC.size());
    glb_sum(CON_F3ygLlgRlgmuRlC_F1xgmuLcC.data(), CON_F3ygLlgRlgmuRlC_F1xgmuLcC.size());
    glb_sum(CON_F3xgmuLlgmuRlgnuLl_F1ygnuLc.data(), CON_F3xgmuLlgmuRlgnuLl_F1ygnuLc.size());
    glb_sum(CON_F3ygnuLlgnuRlgmuLl_F1xgmuLc.data(), CON_F3ygnuLlgnuRlgmuLl_F1xgmuLc.size());
    glb_sum(CON_F3xgmuLlgmuRlgLl_F1ygRl.data(), CON_F3xgmuLlgmuRlgLl_F1ygRl.size());
    glb_sum(CON_F3ygnuLlgnuRlgLl_F1xgRl.data(), CON_F3ygnuLlgnuRlgLl_F1xgRl.size());
    glb_sum(CON_F3xgmuLlgmuRlgRl_F1ygLl.data(), CON_F3xgmuLlgmuRlgRl_F1ygLl.size());
    glb_sum(CON_F3ygnuLlgnuRlgRl_F1xgLl.data(), CON_F3ygnuLlgnuRlgRl_F1xgLl.size());
    glb_sum(CON_F3xgmuLlgmuRlgnuLl_F1ygnuLcC.data(), CON_F3xgmuLlgmuRlgnuLl_F1ygnuLcC.size());
    glb_sum(CON_F3ygnuLlgnuRlgmuLl_F1xgmuLcC.data(), CON_F3ygnuLlgnuRlgmuLl_F1xgmuLcC.size());
    glb_sum(CON_F3xgmuLcgmuLlgRl_F1ygLl.data(), CON_F3xgmuLcgmuLlgRl_F1ygLl.size());
    glb_sum(CON_F3ygnuLcgnuLlgRl_F1xgLl.data(), CON_F3ygnuLcgnuLlgRl_F1xgLl.size());
    glb_sum(CON_F3xgmuRlgmuLlgRl_F1ygLl.data(), CON_F3xgmuRlgmuLlgRl_F1ygLl.size());
    glb_sum(CON_F3ygnuRlgnuLlgRl_F1xgLl.data(), CON_F3ygnuRlgnuLlgRl_F1xgLl.size());
    glb_sum(CON_F3xgLcgRlgRlC_F1ygLl.data(), CON_F3xgLcgRlgRlC_F1ygLl.size());
    glb_sum(CON_F3ygLcgRlgRlC_F1xgLl.data(), CON_F3ygLcgRlgRlC_F1xgLl.size());
    glb_sum(CON_F3xgmuRlgmuLlgnuLl_F1ygnuLc.data(), CON_F3xgmuRlgmuLlgnuLl_F1ygnuLc.size());
    glb_sum(CON_F3ygnuRlgnuLlgmuLl_F1xgmuLc.data(), CON_F3ygnuRlgnuLlgmuLl_F1xgmuLc.size());
    glb_sum(CON_F3xgmuRlgmuLlgLl_F1ygRl.data(), CON_F3xgmuRlgmuLlgLl_F1ygRl.size());
    glb_sum(CON_F3ygnuRlgnuLlgLl_F1xgRl.data(), CON_F3ygnuRlgnuLlgLl_F1xgRl.size());
    glb_sum(CON_F3xgmuRlgmuLlgnuLl_F1ygnuLcC.data(), CON_F3xgmuRlgmuLlgnuLl_F1ygnuLcC.size());
    glb_sum(CON_F3ygnuRlgnuLlgmuLl_F1xgmuLcC.data(), CON_F3ygnuRlgnuLlgmuLl_F1xgmuLcC.size());
    glb_sum(CON_F3xgmuLcgmuLlgLl_F1ygRl.data(), CON_F3xgmuLcgmuLlgLl_F1ygRl.size());
    glb_sum(CON_F3ygnuLcgnuLlgLl_F1xgRl.data(), CON_F3ygnuLcgnuLlgLl_F1xgRl.size());
    glb_sum(CON_F3xgLcgRlgLlC_F1ygRl.data(), CON_F3xgLcgRlgLlC_F1ygRl.size());
    glb_sum(CON_F3ygLcgRlgLlC_F1xgRl.data(), CON_F3ygLcgRlgLlC_F1xgRl.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLl.data(), CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLl.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLc.data(), CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLc.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLl.data(), CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLl.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLlC.data(), CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLlC.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLl.data(), CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLl.size());
    glb_sum(CON_F2gmuLlgLl_F1xgmuLl_F1ygRl.data(), CON_F2gmuLlgLl_F1xgmuLl_F1ygRl.size());
    glb_sum(CON_F2gLlgnuLl_F1xgRl_F1ygnuLl.data(), CON_F2gLlgnuLl_F1xgRl_F1ygnuLl.size());
    glb_sum(CON_F2gmuLlgRl_F1xgmuLl_F1ygLl.data(), CON_F2gmuLlgRl_F1xgmuLl_F1ygLl.size());
    glb_sum(CON_F2gRlgnuLl_F1xgLl_F1ygnuLl.data(), CON_F2gRlgnuLl_F1xgLl_F1ygnuLl.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLcC.data(), CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLcC.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLl.data(), CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLl.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLc.data(), CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLc.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLlC.data(), CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLlC.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLlC.data(), CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLlC.size());
    glb_sum(CON_F2gmuLlgLl_F1xgmuLlC_F1ygRl.data(), CON_F2gmuLlgLl_F1xgmuLlC_F1ygRl.size());
    glb_sum(CON_F2gLlgnuLl_F1xgRl_F1ygnuLlC.data(), CON_F2gLlgnuLl_F1xgRl_F1ygnuLlC.size());
    glb_sum(CON_F2gmuLlgRl_F1xgmuLlC_F1ygLl.data(), CON_F2gmuLlgRl_F1xgmuLlC_F1ygLl.size());
    glb_sum(CON_F2gRlgnuLl_F1xgLl_F1ygnuLlC.data(), CON_F2gRlgnuLl_F1xgLl_F1ygnuLlC.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLcC.data(), CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLcC.size());
    glb_sum(CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLlC.data(), CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLlC.size());
    glb_sum(CON_F2gmuLlgRl_F1xgmuLc_F1ygLl.data(), CON_F2gmuLlgRl_F1xgmuLc_F1ygLl.size());
    glb_sum(CON_F2gRlgnuLl_F1xgLl_F1ygnuLc.data(), CON_F2gRlgnuLl_F1xgLl_F1ygnuLc.size());
    glb_sum(CON_F2gLlgRl_F1xgRl_F1ygLl.data(), CON_F2gLlgRl_F1xgRl_F1ygLl.size());
    glb_sum(CON_F2gRlgLl_F1xgLl_F1ygRl.data(), CON_F2gRlgLl_F1xgLl_F1ygRl.size());
    glb_sum(CON_F2gRlgRl_F1xgLl_F1ygLl.data(), CON_F2gRlgRl_F1xgLl_F1ygLl.size());
    glb_sum(CON_F2gmuLlgRl_F1xgmuLcC_F1ygLl.data(), CON_F2gmuLlgRl_F1xgmuLcC_F1ygLl.size());
    glb_sum(CON_F2gRlgnuLl_F1xgLl_F1ygnuLcC.data(), CON_F2gRlgnuLl_F1xgLl_F1ygnuLcC.size());
    glb_sum(CON_F2gmuLlgLl_F1xgmuLc_F1ygRl.data(), CON_F2gmuLlgLl_F1xgmuLc_F1ygRl.size());
    glb_sum(CON_F2gLlgnuLl_F1xgRl_F1ygnuLc.data(), CON_F2gLlgnuLl_F1xgRl_F1ygnuLc.size());
    glb_sum(CON_F2gLlgLl_F1xgRl_F1ygRl.data(), CON_F2gLlgLl_F1xgRl_F1ygRl.size());
    glb_sum(CON_F2gmuLlgLl_F1xgmuLcC_F1ygRl.data(), CON_F2gmuLlgLl_F1xgmuLcC_F1ygRl.size());
    glb_sum(CON_F2gLlgnuLl_F1xgRl_F1ygnuLcC.data(), CON_F2gLlgnuLl_F1xgRl_F1ygnuLcC.size());
    glb_sum(CON_F3xgmuLlgmuLlgRl.data(), CON_F3xgmuLlgmuLlgRl.size());
    glb_sum(CON_F3ygnuLlgnuLlgRl.data(), CON_F3ygnuLlgnuLlgRl.size());
    glb_sum(CON_F3xgmuLlgmuLlgLl.data(), CON_F3xgmuLlgmuLlgLl.size());
    glb_sum(CON_F3ygnuLlgnuLlgLl.data(), CON_F3ygnuLlgnuLlgLl.size());
    glb_sum(CON_F3xgLlgRlgRl.data(), CON_F3xgLlgRlgRl.size());
    glb_sum(CON_F3ygLlgRlgRl.data(), CON_F3ygLlgRlgRl.size());
    glb_sum(CON_F3xgLlgRlgLl.data(), CON_F3xgLlgRlgLl.size());
    glb_sum(CON_F3ygLlgRlgLl.data(), CON_F3ygLlgRlgLl.size());
    glb_sum(CON_F3xgmuLlgmuRlgRl.data(), CON_F3xgmuLlgmuRlgRl.size());
    glb_sum(CON_F3ygnuLlgnuRlgRl.data(), CON_F3ygnuLlgnuRlgRl.size());
    glb_sum(CON_F3xgmuLlgmuRlgLl.data(), CON_F3xgmuLlgmuRlgLl.size());
    glb_sum(CON_F3ygnuLlgnuRlgLl.data(), CON_F3ygnuLlgnuRlgLl.size());
    glb_sum(CON_F3xgmuRlgmuLlgRl.data(), CON_F3xgmuRlgmuLlgRl.size());
    glb_sum(CON_F3ygnuRlgnuLlgRl.data(), CON_F3ygnuRlgnuLlgRl.size());
    glb_sum(CON_F3xgmuRlgmuLlgLl.data(), CON_F3xgmuRlgmuLlgLl.size());
    glb_sum(CON_F3ygnuRlgnuLlgLl.data(), CON_F3ygnuRlgnuLlgLl.size());
    glb_sum(CON_F3xgmuLcgmuLlgLl.data(), CON_F3xgmuLcgmuLlgLl.size());
    glb_sum(CON_F3ygnuLcgnuLlgLl.data(), CON_F3ygnuLcgnuLlgLl.size());
    glb_sum(CON_F3xgLcgRlgLl.data(), CON_F3xgLcgRlgLl.size());
    glb_sum(CON_F3ygLcgRlgLl.data(), CON_F3ygLcgRlgLl.size());
    glb_sum(CON_F3xgmuLcgmuLlgRl.data(), CON_F3xgmuLcgmuLlgRl.size());
    glb_sum(CON_F3ygnuLcgnuLlgRl.data(), CON_F3ygnuLcgnuLlgRl.size());
    glb_sum(CON_F3xgLcgRlgRl.data(), CON_F3xgLcgRlgRl.size());
    glb_sum(CON_F3ygLcgRlgRl.data(), CON_F3ygLcgRlgRl.size());
    glb_sum(CON_F2gmuLlgRl_F1xgmuLl.data(), CON_F2gmuLlgRl_F1xgmuLl.size());
    glb_sum(CON_F2gRlgnuLl_F1ygnuLl.data(), CON_F2gRlgnuLl_F1ygnuLl.size());
    glb_sum(CON_F2gmuLlgLl_F1xgmuLl.data(), CON_F2gmuLlgLl_F1xgmuLl.size());
    glb_sum(CON_F2gLlgnuLl_F1ygnuLl.data(), CON_F2gLlgnuLl_F1ygnuLl.size());
    glb_sum(CON_F2gmuLlgRl_F1xgmuLlC.data(), CON_F2gmuLlgRl_F1xgmuLlC.size());
    glb_sum(CON_F2gRlgnuLl_F1ygnuLlC.data(), CON_F2gRlgnuLl_F1ygnuLlC.size());
    glb_sum(CON_F2gmuLlgLl_F1xgmuLlC.data(), CON_F2gmuLlgLl_F1xgmuLlC.size());
    glb_sum(CON_F2gLlgnuLl_F1ygnuLlC.data(), CON_F2gLlgnuLl_F1ygnuLlC.size());
    glb_sum(CON_F2gRlgRl_F1xgLl.data(), CON_F2gRlgRl_F1xgLl.size());
    glb_sum(CON_F2gRlgRl_F1ygLl.data(), CON_F2gRlgRl_F1ygLl.size());
    glb_sum(CON_F2gRlgLl_F1xgLl.data(), CON_F2gRlgLl_F1xgLl.size());
    glb_sum(CON_F2gLlgRl_F1ygLl.data(), CON_F2gLlgRl_F1ygLl.size());
    glb_sum(CON_F2gLlgRl_F1xgRl.data(), CON_F2gLlgRl_F1xgRl.size());
    glb_sum(CON_F2gRlgLl_F1ygRl.data(), CON_F2gRlgLl_F1ygRl.size());
    glb_sum(CON_F2gLlgLl_F1xgRl.data(), CON_F2gLlgLl_F1xgRl.size());
    glb_sum(CON_F2gLlgLl_F1ygRl.data(), CON_F2gLlgLl_F1ygRl.size());
    glb_sum(CON_F2gmuLlgLl_F1xgmuLc.data(), CON_F2gmuLlgLl_F1xgmuLc.size());
    glb_sum(CON_F2gLlgnuLl_F1ygnuLc.data(), CON_F2gLlgnuLl_F1ygnuLc.size());
    glb_sum(CON_F2gmuLlgLl_F1xgmuLcC.data(), CON_F2gmuLlgLl_F1xgmuLcC.size());
    glb_sum(CON_F2gLlgnuLl_F1ygnuLcC.data(), CON_F2gLlgnuLl_F1ygnuLcC.size());
    glb_sum(CON_F2gmuLlgRl_F1xgmuLc.data(), CON_F2gmuLlgRl_F1xgmuLc.size());
    glb_sum(CON_F2gRlgnuLl_F1ygnuLc.data(), CON_F2gRlgnuLl_F1ygnuLc.size());
    glb_sum(CON_F2gmuLlgRl_F1xgmuLcC.data(), CON_F2gmuLlgRl_F1xgmuLcC.size());
    glb_sum(CON_F2gRlgnuLl_F1ygnuLcC.data(), CON_F2gRlgnuLl_F1ygnuLcC.size());
    glb_sum(CON_F2gLlgRl.data(), CON_F2gLlgRl.size());
    glb_sum(CON_F2gRlgLl.data(), CON_F2gRlgLl.size());
    glb_sum(CON_F2gLlgLl.data(), CON_F2gLlgLl.size());
    glb_sum(CON_F2gRlgRl.data(), CON_F2gRlgRl.size());
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
    if ( 0 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gmuLlgnuLlgmuLlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gmuLlgnuLlgmuLlgnuLl[pt] /= nsrc;
        auto temp = CON_F4gmuLlgnuLlgmuLlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 1 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gmuLlgLlgmuLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gmuLlgLlgmuLlgRl[pt] /= nsrc;
        auto temp = CON_F4gmuLlgLlgmuLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 2 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gLlgnuLlgRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gLlgnuLlgRlgnuLl[pt] /= nsrc;
        auto temp = CON_F4gLlgnuLlgRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 3 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gmuLlgnuLlgmuLlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gmuLlgnuLlgmuLlgnuRl[pt] /= nsrc;
        auto temp = CON_F4gmuLlgnuLlgmuLlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 4 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gmuLlgnuLlgmuRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gmuLlgnuLlgmuRlgnuLl[pt] /= nsrc;
        auto temp = CON_F4gmuLlgnuLlgmuRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 5 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gmuLlgRlgmuRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gmuLlgRlgmuRlgLl[pt] /= nsrc;
        auto temp = CON_F4gmuLlgRlgmuRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 6 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gLlgnuLlgRlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gLlgnuLlgRlgnuRl[pt] /= nsrc;
        auto temp = CON_F4gLlgnuLlgRlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 7 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gmuLlgnuLlgmuRlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gmuLlgnuLlgmuRlgnuRl[pt] /= nsrc;
        auto temp = CON_F4gmuLlgnuLlgmuRlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 8 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4gLlgLlgRlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4gLlgLlgRlgRl[pt] /= nsrc;
        auto temp = CON_F4gLlgLlgRlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 9 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F2gmuLlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F2gmuLlgnuLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F2gmuLlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 10 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F2gmuLlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F2gmuLlgnuRl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F2gmuLlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 11 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F2gmuLlgnuRlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F2gmuLlgnuRlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F2gmuLlgnuRlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 12 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F2gmuLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F2gmuLlgRl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F2gmuLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 13 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F2gRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F2gRlgnuLl[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F2gRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 14 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F2gmuLlgnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F2gmuLlgnuLlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F2gmuLlgnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 15 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgRl_F2gRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgRl_F2gRlgLl[pt] /= nsrc;
        auto temp = CON_F2gLlgRl_F2gRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 16 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F2gmuLlgRlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F2gmuLlgRlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F2gmuLlgRlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 17 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F2gRlgnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F2gRlgnuLlC[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F2gRlgnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 18 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuLlgnuLlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuLlgnuLlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuLlgnuLlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 19 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuLlgnuLcgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuLlgnuLcgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuLlgnuLcgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 20 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLcgmuLlgnuLlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLcgmuLlgnuLlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLcgmuLlgnuLlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 21 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuLlgRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuLlgRlgLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuLlgRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 22 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLlgRlgnuRlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLlgRlgnuRlgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgLlgRlgnuRlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 23 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuLlgnuRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuLlgnuRlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuLlgnuRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 24 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuRlgmuLlgnuLlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuRlgmuLlgnuLlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuRlgmuLlgnuLlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 25 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuLlgnuLlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuLlgnuLlgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuLlgnuLlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 26 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuRlgnuLlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuRlgnuLlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuRlgnuLlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 27 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuLlgRcgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuLlgRcgLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuLlgRcgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 28 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLcgRlgnuRlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLcgRlgnuRlgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgLcgRlgnuRlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 29 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLlgRlgnuRcgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLlgRlgnuRcgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgLlgRlgnuRcgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 30 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLcgmuLlgRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLcgmuLlgRlgLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLcgmuLlgRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 31 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLlgRlgLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLlgRlgLlgRl[pt] /= nsrc;
        auto temp = CON_F4pgLlgRlgLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 32 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuRlgmuLlgLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuRlgmuLlgLlgRl[pt] /= nsrc;
        auto temp = CON_F4pgmuRlgmuLlgLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 33 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLlgRlgnuRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLlgRlgnuRlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgLlgRlgnuRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 34 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuRlgLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuRlgLlgRl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuRlgLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 35 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLlgRlgnuLlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLlgRlgnuLlgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgLlgRlgnuLlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 36 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLlgRlgLcgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLlgRlgLcgRl[pt] /= nsrc;
        auto temp = CON_F4pgLlgRlgLcgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 37 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLcgRlgLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLcgRlgLlgRl[pt] /= nsrc;
        auto temp = CON_F4pgLcgRlgLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 38 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuRlgnuLcgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuRlgnuLcgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuRlgnuLcgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 39 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLcgmuLlgnuLlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLcgmuLlgnuLlgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgmuLcgmuLlgnuLlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 40 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuRlgnuRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuRlgnuRlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuRlgnuRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 41 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuRlgmuLlgnuLlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuRlgmuLlgnuLlgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgmuRlgmuLlgnuLlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 42 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuRlgnuLlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuRlgnuLlgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuRlgnuLlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 43 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLlgmuRlgLcgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLlgmuRlgLcgRl[pt] /= nsrc;
        auto temp = CON_F4pgmuLlgmuRlgLcgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 44 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLcgRlgnuLlgnuRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLcgRlgnuLlgnuRl[pt] /= nsrc;
        auto temp = CON_F4pgLcgRlgnuLlgnuRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 45 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuRlgmuLlgnuLcgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuRlgmuLlgnuLcgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuRlgmuLlgnuLcgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 46 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuLcgmuLlgnuRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuLcgmuLlgnuRlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuLcgmuLlgnuRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 47 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuRlgmuLlgnuRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuRlgmuLlgnuRlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgmuRlgmuLlgnuRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 48 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgmuRlgmuLlgLcgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgmuRlgmuLlgLcgRl[pt] /= nsrc;
        auto temp = CON_F4pgmuRlgmuLlgLcgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 49 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F4pgLcgRlgnuRlgnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F4pgLcgRlgnuRlgnuLl[pt] /= nsrc;
        auto temp = CON_F4pgLcgRlgnuRlgnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 50 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuLlgnuLl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuLlgnuLl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuLlgnuLl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 51 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuLlgmuLl_F1xgmuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuLlgmuLl_F1xgmuLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuLlgmuLl_F1xgmuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 52 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLcgmuLlgnuLl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLcgmuLlgnuLl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLcgmuLlgnuLl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 53 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLcgnuLlgmuLl_F1xgmuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLcgnuLlgmuLl_F1xgmuLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLcgnuLlgmuLl_F1xgmuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 54 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLlgRlgnuRlC_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLlgRlgnuRlC_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F3xgLlgRlgnuRlC_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 55 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLlgRlgmuRlC_F1xgmuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLlgRlgmuRlC_F1xgmuLl[pt] /= nsrc;
        auto temp = CON_F3ygLlgRlgmuRlC_F1xgmuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 56 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuRlgmuLlgnuLl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuRlgmuLlgnuLl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F3xgmuRlgmuLlgnuLl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 57 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuRlgnuLlgmuLl_F1xgmuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuRlgnuLlgmuLl_F1xgmuLl[pt] /= nsrc;
        auto temp = CON_F3ygnuRlgnuLlgmuLl_F1xgmuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 58 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuRlgnuLl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuRlgnuLl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuRlgnuLl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 59 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuRlgmuLl_F1xgmuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuRlgmuLl_F1xgmuLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuRlgmuLl_F1xgmuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 60 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLcgRlgnuRlC_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLcgRlgnuRlC_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F3xgLcgRlgnuRlC_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 61 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLcgRlgmuRlC_F1xgmuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLcgRlgmuRlC_F1xgmuLl[pt] /= nsrc;
        auto temp = CON_F3ygLcgRlgmuRlC_F1xgmuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 62 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuLlgnuLl_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuLlgnuLl_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuLlgnuLl_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 63 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuLlgmuLl_F1xgmuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuLlgmuLl_F1xgmuLc[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuLlgmuLl_F1xgmuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 64 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuLlgnuLl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuLlgnuLl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuLlgnuLl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 65 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuLlgmuLl_F1xgmuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuLlgmuLl_F1xgmuLlC[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuLlgmuLl_F1xgmuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 66 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuLlgLl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuLlgLl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuLlgLl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 67 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuLlgLl_F1xgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuLlgLl_F1xgRl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuLlgLl_F1xgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 68 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuLlgRl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuLlgRl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuLlgRl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 69 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuLlgRl_F1xgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuLlgRl_F1xgLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuLlgRl_F1xgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 70 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuLlgnuLl_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuLlgnuLl_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuLlgnuLl_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 71 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuLlgmuLl_F1xgmuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuLlgmuLl_F1xgmuLcC[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuLlgmuLl_F1xgmuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 72 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLcgmuLlgnuLl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLcgmuLlgnuLl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F3xgmuLcgmuLlgnuLl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 73 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLcgnuLlgmuLl_F1xgmuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLcgnuLlgmuLl_F1xgmuLlC[pt] /= nsrc;
        auto temp = CON_F3ygnuLcgnuLlgmuLl_F1xgmuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 74 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLlgRlgnuRlC_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLlgRlgnuRlC_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F3xgLlgRlgnuRlC_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 75 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLlgRlgmuRlC_F1xgmuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLlgRlgmuRlC_F1xgmuLlC[pt] /= nsrc;
        auto temp = CON_F3ygLlgRlgmuRlC_F1xgmuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 76 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuRlgmuLlgnuLl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuRlgmuLlgnuLl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F3xgmuRlgmuLlgnuLl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 77 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuRlgnuLlgmuLl_F1xgmuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuRlgnuLlgmuLl_F1xgmuLlC[pt] /= nsrc;
        auto temp = CON_F3ygnuRlgnuLlgmuLl_F1xgmuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 78 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuRlgnuLl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuRlgnuLl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuRlgnuLl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 79 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuRlgmuLl_F1xgmuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuRlgmuLl_F1xgmuLlC[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuRlgmuLl_F1xgmuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 80 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLcgRlgnuRlC_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLcgRlgnuRlC_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F3xgLcgRlgnuRlC_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 81 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLcgRlgmuRlC_F1xgmuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLcgRlgmuRlC_F1xgmuLlC[pt] /= nsrc;
        auto temp = CON_F3ygLcgRlgmuRlC_F1xgmuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 82 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLlgRlgnuRlC_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLlgRlgnuRlC_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F3xgLlgRlgnuRlC_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 83 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLlgRlgmuRlC_F1xgmuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLlgRlgmuRlC_F1xgmuLc[pt] /= nsrc;
        auto temp = CON_F3ygLlgRlgmuRlC_F1xgmuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 84 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLlgRlgLlC_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLlgRlgLlC_F1ygRl[pt] /= nsrc;
        auto temp = CON_F3xgLlgRlgLlC_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 85 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLlgRlgLlC_F1xgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLlgRlgLlC_F1xgRl[pt] /= nsrc;
        auto temp = CON_F3ygLlgRlgLlC_F1xgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 86 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLlgRlgRlC_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLlgRlgRlC_F1ygLl[pt] /= nsrc;
        auto temp = CON_F3xgLlgRlgRlC_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 87 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLlgRlgRlC_F1xgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLlgRlgRlC_F1xgLl[pt] /= nsrc;
        auto temp = CON_F3ygLlgRlgRlC_F1xgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 88 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLlgRlgnuRlC_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLlgRlgnuRlC_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F3xgLlgRlgnuRlC_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 89 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLlgRlgmuRlC_F1xgmuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLlgRlgmuRlC_F1xgmuLcC[pt] /= nsrc;
        auto temp = CON_F3ygLlgRlgmuRlC_F1xgmuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 90 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuRlgnuLl_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuRlgnuLl_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuRlgnuLl_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 91 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuRlgmuLl_F1xgmuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuRlgmuLl_F1xgmuLc[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuRlgmuLl_F1xgmuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 92 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuRlgLl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuRlgLl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuRlgLl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 93 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuRlgLl_F1xgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuRlgLl_F1xgRl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuRlgLl_F1xgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 94 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuRlgRl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuRlgRl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuRlgRl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 95 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuRlgRl_F1xgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuRlgRl_F1xgLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuRlgRl_F1xgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 96 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuRlgnuLl_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuRlgnuLl_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuRlgnuLl_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 97 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuRlgmuLl_F1xgmuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuRlgmuLl_F1xgmuLcC[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuRlgmuLl_F1xgmuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 98 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLcgmuLlgRl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLcgmuLlgRl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLcgmuLlgRl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 99 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLcgnuLlgRl_F1xgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLcgnuLlgRl_F1xgLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLcgnuLlgRl_F1xgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 100 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuRlgmuLlgRl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuRlgmuLlgRl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F3xgmuRlgmuLlgRl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 101 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuRlgnuLlgRl_F1xgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuRlgnuLlgRl_F1xgLl[pt] /= nsrc;
        auto temp = CON_F3ygnuRlgnuLlgRl_F1xgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 102 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLcgRlgRlC_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLcgRlgRlC_F1ygLl[pt] /= nsrc;
        auto temp = CON_F3xgLcgRlgRlC_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 103 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLcgRlgRlC_F1xgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLcgRlgRlC_F1xgLl[pt] /= nsrc;
        auto temp = CON_F3ygLcgRlgRlC_F1xgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 104 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuRlgmuLlgnuLl_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuRlgmuLlgnuLl_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F3xgmuRlgmuLlgnuLl_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 105 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuRlgnuLlgmuLl_F1xgmuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuRlgnuLlgmuLl_F1xgmuLc[pt] /= nsrc;
        auto temp = CON_F3ygnuRlgnuLlgmuLl_F1xgmuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 106 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuRlgmuLlgLl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuRlgmuLlgLl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F3xgmuRlgmuLlgLl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 107 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuRlgnuLlgLl_F1xgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuRlgnuLlgLl_F1xgRl[pt] /= nsrc;
        auto temp = CON_F3ygnuRlgnuLlgLl_F1xgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 108 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuRlgmuLlgnuLl_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuRlgmuLlgnuLl_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F3xgmuRlgmuLlgnuLl_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 109 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuRlgnuLlgmuLl_F1xgmuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuRlgnuLlgmuLl_F1xgmuLcC[pt] /= nsrc;
        auto temp = CON_F3ygnuRlgnuLlgmuLl_F1xgmuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 110 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLcgmuLlgLl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLcgmuLlgLl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F3xgmuLcgmuLlgLl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 111 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLcgnuLlgLl_F1xgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLcgnuLlgLl_F1xgRl[pt] /= nsrc;
        auto temp = CON_F3ygnuLcgnuLlgLl_F1xgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 112 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLcgRlgLlC_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLcgRlgLlC_F1ygRl[pt] /= nsrc;
        auto temp = CON_F3xgLcgRlgLlC_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 113 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLcgRlgLlC_F1xgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLcgRlgLlC_F1xgRl[pt] /= nsrc;
        auto temp = CON_F3ygLcgRlgLlC_F1xgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 114 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 115 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 116 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 117 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 118 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 119 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F1xgmuLl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F1xgmuLl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F1xgmuLl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 120 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F1xgRl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F1xgRl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F1xgRl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 121 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgRl_F1xgmuLl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgRl_F1xgmuLl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgRl_F1xgmuLl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 122 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgnuLl_F1xgLl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgnuLl_F1xgLl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F2gRlgnuLl_F1xgLl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 123 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLl_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 124 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 125 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 126 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLc_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 127 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 128 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F1xgmuLlC_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F1xgmuLlC_F1ygRl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F1xgmuLlC_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 129 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F1xgRl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F1xgRl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F1xgRl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 130 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgRl_F1xgmuLlC_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgRl_F1xgmuLlC_F1ygLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgRl_F1xgmuLlC_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 131 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgnuLl_F1xgLl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgnuLl_F1xgLl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F2gRlgnuLl_F1xgLl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 132 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLlC_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 133 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgnuLl_F1xgmuLcC_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 134 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgRl_F1xgmuLc_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgRl_F1xgmuLc_F1ygLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgRl_F1xgmuLc_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 135 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgnuLl_F1xgLl_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgnuLl_F1xgLl_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F2gRlgnuLl_F1xgLl_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 136 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgRl_F1xgRl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgRl_F1xgRl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F2gLlgRl_F1xgRl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 137 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgLl_F1xgLl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgLl_F1xgLl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F2gRlgLl_F1xgLl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 138 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgRl_F1xgLl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgRl_F1xgLl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F2gRlgRl_F1xgLl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 139 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgRl_F1xgmuLcC_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgRl_F1xgmuLcC_F1ygLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgRl_F1xgmuLcC_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 140 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgnuLl_F1xgLl_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgnuLl_F1xgLl_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F2gRlgnuLl_F1xgLl_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 141 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F1xgmuLc_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F1xgmuLc_F1ygRl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F1xgmuLc_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 142 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F1xgRl_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F1xgRl_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F1xgRl_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 143 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgLl_F1xgRl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgLl_F1xgRl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F2gLlgLl_F1xgRl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 144 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F1xgmuLcC_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F1xgmuLcC_F1ygRl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F1xgmuLcC_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 145 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F1xgRl_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F1xgRl_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F1xgRl_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 146 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuLlgRl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 147 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuLlgRl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 148 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuLlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuLlgLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuLlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 149 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuLlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuLlgLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuLlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 150 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLlgRlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLlgRlgRl[pt] /= nsrc;
        auto temp = CON_F3xgLlgRlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 151 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLlgRlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLlgRlgRl[pt] /= nsrc;
        auto temp = CON_F3ygLlgRlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 152 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLlgRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLlgRlgLl[pt] /= nsrc;
        auto temp = CON_F3xgLlgRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 153 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLlgRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLlgRlgLl[pt] /= nsrc;
        auto temp = CON_F3ygLlgRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 154 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuRlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuRlgRl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuRlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 155 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuRlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuRlgRl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuRlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 156 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLlgmuRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLlgmuRlgLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLlgmuRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 157 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLlgnuRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLlgnuRlgLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLlgnuRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 158 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuRlgmuLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuRlgmuLlgRl[pt] /= nsrc;
        auto temp = CON_F3xgmuRlgmuLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 159 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuRlgnuLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuRlgnuLlgRl[pt] /= nsrc;
        auto temp = CON_F3ygnuRlgnuLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 160 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuRlgmuLlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuRlgmuLlgLl[pt] /= nsrc;
        auto temp = CON_F3xgmuRlgmuLlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 161 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuRlgnuLlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuRlgnuLlgLl[pt] /= nsrc;
        auto temp = CON_F3ygnuRlgnuLlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 162 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLcgmuLlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLcgmuLlgLl[pt] /= nsrc;
        auto temp = CON_F3xgmuLcgmuLlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 163 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLcgnuLlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLcgnuLlgLl[pt] /= nsrc;
        auto temp = CON_F3ygnuLcgnuLlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 164 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLcgRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLcgRlgLl[pt] /= nsrc;
        auto temp = CON_F3xgLcgRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 165 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLcgRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLcgRlgLl[pt] /= nsrc;
        auto temp = CON_F3ygLcgRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 166 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgmuLcgmuLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgmuLcgmuLlgRl[pt] /= nsrc;
        auto temp = CON_F3xgmuLcgmuLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 167 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygnuLcgnuLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygnuLcgnuLlgRl[pt] /= nsrc;
        auto temp = CON_F3ygnuLcgnuLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 168 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3xgLcgRlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3xgLcgRlgRl[pt] /= nsrc;
        auto temp = CON_F3xgLcgRlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 169 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F3ygLcgRlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F3ygLcgRlgRl[pt] /= nsrc;
        auto temp = CON_F3ygLcgRlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 170 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgRl_F1xgmuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgRl_F1xgmuLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgRl_F1xgmuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 171 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgnuLl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgnuLl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F2gRlgnuLl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 172 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F1xgmuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F1xgmuLl[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F1xgmuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 173 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F1ygnuLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F1ygnuLl[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F1ygnuLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 174 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgRl_F1xgmuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgRl_F1xgmuLlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgRl_F1xgmuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 175 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgnuLl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgnuLl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F2gRlgnuLl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 176 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F1xgmuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F1xgmuLlC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F1xgmuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 177 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F1ygnuLlC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F1ygnuLlC[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F1ygnuLlC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 178 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgRl_F1xgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgRl_F1xgLl[pt] /= nsrc;
        auto temp = CON_F2gRlgRl_F1xgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 179 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgRl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgRl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F2gRlgRl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 180 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgLl_F1xgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgLl_F1xgLl[pt] /= nsrc;
        auto temp = CON_F2gRlgLl_F1xgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 181 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgRl_F1ygLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgRl_F1ygLl[pt] /= nsrc;
        auto temp = CON_F2gLlgRl_F1ygLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 182 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgRl_F1xgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgRl_F1xgRl[pt] /= nsrc;
        auto temp = CON_F2gLlgRl_F1xgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 183 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgLl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgLl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F2gRlgLl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 184 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgLl_F1xgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgLl_F1xgRl[pt] /= nsrc;
        auto temp = CON_F2gLlgLl_F1xgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 185 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgLl_F1ygRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgLl_F1ygRl[pt] /= nsrc;
        auto temp = CON_F2gLlgLl_F1ygRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 186 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F1xgmuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F1xgmuLc[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F1xgmuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 187 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 188 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgLl_F1xgmuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgLl_F1xgmuLcC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgLl_F1xgmuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 189 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgnuLl_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgnuLl_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F2gLlgnuLl_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 190 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgRl_F1xgmuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgRl_F1xgmuLc[pt] /= nsrc;
        auto temp = CON_F2gmuLlgRl_F1xgmuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 191 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgnuLl_F1ygnuLc", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgnuLl_F1ygnuLc[pt] /= nsrc;
        auto temp = CON_F2gRlgnuLl_F1ygnuLc[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 192 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gmuLlgRl_F1xgmuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gmuLlgRl_F1xgmuLcC[pt] /= nsrc;
        auto temp = CON_F2gmuLlgRl_F1xgmuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 193 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgnuLl_F1ygnuLcC", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgnuLl_F1ygnuLcC[pt] /= nsrc;
        auto temp = CON_F2gRlgnuLl_F1ygnuLcC[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 194 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgRl[pt] /= nsrc;
        auto temp = CON_F2gLlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 195 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgLl[pt] /= nsrc;
        auto temp = CON_F2gRlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 196 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gLlgLl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gLlgLl[pt] /= nsrc;
        auto temp = CON_F2gLlgLl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }
    if ( 197 % nnode == nodeid ) {
      char filename[1024];
      sprintf(filename, "Output/traj_%d_F2gRlgRl", meas_arg.TrajCur);
      file.open(filename,std::ios::out | std::ios::binary);
      for ( long long int pt = 0; pt < glb_v; pt++ ) {
        CON_F2gRlgRl[pt] /= nsrc;
        auto temp = CON_F2gRlgRl[pt].real();
        file.write((const char*)&temp, 8);
      }
      file.close();
      printf ("Binary file saved on %s\n",filename);
    }

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



