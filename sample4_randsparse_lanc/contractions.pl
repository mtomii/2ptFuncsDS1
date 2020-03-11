#!/usr/bin/perl
use strict;
use warnings;
use Functions;
use CorrManage;
use FormatManage;

my $indiv_flag = 0;# 1 for checking consistency up to random numbers of multiple point sources

my $tmp_meas_code = "./main.tmp.C";
my $meas_code = "./main.C";
my $nr = 20;
my $nsrc_s = 4;
my $nsrc_t = 4;

my $ch_list = "channels.list";
my $decomp_dir = "./decompositions/";

my $acls_type_4 = &FormatManage::get_type_4();
my $acls_type_22 = &FormatManage::get_type_22();
my $acls_type_4p = &FormatManage::get_type_4p();
my $acls_type_31 = &FormatManage::get_type_31();
my $acls_type_211 = &FormatManage::get_type_211();
my $acls_type_3 = &FormatManage::get_type_3();
my $acls_type_21 = &FormatManage::get_type_21();
my $acls_type_2 = &FormatManage::get_type_2();
my @acls_types = ($acls_type_4,$acls_type_22,$acls_type_4p,
                  $acls_type_31,$acls_type_211,$acls_type_3,
                  $acls_type_21,$acls_type_2);

my $rtype_4    = &RepManage::get_type_4();
my $rtype_22   = &RepManage::get_type_22();
my $rtype_4p   = &RepManage::get_type_4p();
my $rtype_31xy = &RepManage::get_type_31xy();
my $rtype_31yx = &RepManage::get_type_31yx();
my $rtype_211  = &RepManage::get_type_211();
my $rtype_3x   = &RepManage::get_type_3x();
my $rtype_3y   = &RepManage::get_type_3y();
my $rtype_21x  = &RepManage::get_type_21x();
my $rtype_21y  = &RepManage::get_type_21y();
my $rtype_2    = &RepManage::get_type_2();

my $tr_type_4  = &TraceManage::get_type_4();
my $tr_type_4p = &TraceManage::get_type_4p();
my $tr_type_3x = &TraceManage::get_type_3x();
my $tr_type_3y = &TraceManage::get_type_3y();
my $tr_type_2  = &TraceManage::get_type_2();
my $tr_type_1x = &TraceManage::get_type_1x();
my $tr_type_1y = &TraceManage::get_type_1y();
my @tr_types = &TraceManage::get_types();

my @tr_ttypes = &TraceManage::get_ttypes();

my $uuL   = " sb_1^1 gmu1L^1^2 dv_1^2 ub_2^3 gmu1L^3^4 uv_2^4";
my $uuLt  = " sb_1^1 gmu1L^1^2 uv_1^2 ub_2^3 gmu1L^3^4 dv_2^4";
my $ddL   = " sb_1^1 gmu1L^1^2 dv_1^2 db_2^3 gmu1L^3^4 dv_2^4";
my $ssL   = " sb_1^1 gmu1L^1^2 dv_1^2 sb_2^3 gmu1L^3^4 sv_2^4";
my $ccL   = " sb_1^1 gmu1L^1^2 dv_1^2 cb_2^3 gmu1L^3^4 cv_2^4";
my $ccLt  = " sb_1^1 gmu1L^1^2 cv_1^2 cb_2^3 gmu1L^3^4 dv_2^4";
my $uuR   = " sb_1^1 gmu1L^1^2 dv_1^2 ub_2^3 gmu1R^3^4 uv_2^4";
my $uuRt  = " - 2 sb_1^1 gR^1^2 uv_1^2 ub_2^3 gL^3^4 dv_2^4";
my $ddR   = " sb_1^1 gmu1L^1^2 dv_1^2 db_2^3 gmu1R^3^4 dv_2^4";
my $ddRt  = " - 2 sb_1^1 gR^1^2 dv_1^2 db_2^3 gL^3^4 dv_2^4";
my $ssR   = " sb_1^1 gmu1L^1^2 dv_1^2 sb_2^3 gmu1R^3^4 sv_2^4";
my $ssRt  = " - 2 sb_1^1 gR^1^2 sv_1^2 sb_2^3 gL^3^4 dv_2^4";
my $ccR   = " sb_1^1 gmu1L^1^2 dv_1^2 cb_2^3 gmu1R^3^4 cv_2^4";
my $ccRt  = " - 2 sb_1^1 gR^1^2 cv_1^2 cb_2^3 gL^3^4 dv_2^4";
my $sdL   = " sb_1^1 gL^1^2 dv_1^2";
my $sdR   = " sb_1^1 gR^1^2 dv_1^2";

# F(x)
#my @Fx =  ($uuL,$uuLt,$ddL,      $ssL,      $ccL,$ccLt,
#           $uuR,$uuRt,$ddR,$ddRt,$ssR,$ssRt,$ccR,$ccRt);
my @Fx =  ($uuL,$uuLt,$ddL,      $ssL,
           $uuR,$uuRt,$ddR,$ddRt,$ssR,$ssRt,
           $sdL,$sdR);
my @nameFx = ("uuL","uuLm","ddL",       "ssL",
              "uuR","uuRm","ddR","ddRm","ssR","ssRm",
              "blL","blR");

my $OuuL   = "[(\\bar s_\\alpha d_\\alpha)_L(\\bar u_\\beta u_\\beta)_L](x)";
my $OuuLt  = "[(\\bar s_\\alpha d_\\beta)_L(\\bar u_\\beta u_\\alpha)_L](x)";
my $OddL   = "[(\\bar s_\\alpha d_\\alpha)_L(\\bar d_\\beta d_\\beta)_L](x)";
my $OssL   = "[(\\bar s_\\alpha d_\\alpha)_L(\\bar s_\\beta s_\\beta)_L](x)";
my $OccL   = "[(\\bar s_\\alpha d_\\alpha)_L(\\bar c_\\beta c_\\beta)_L](x)";
my $OccLt  = "[(\\bar s_\\alpha d_\\beta)_L(\\bar c_\\beta c_\\alpha)_L](x)";
my $OuuR   = "[(\\bar s_\\alpha d_\\alpha)_L(\\bar u_\\beta u_\\beta)_R](x)";
my $OuuRt  = "[(\\bar s_\\alpha d_\\beta)_L(\\bar u_\\beta u_\\alpha)_R](x)";
my $OddR   = "[(\\bar s_\\alpha d_\\alpha)_L(\\bar d_\\beta d_\\beta)_R](x)";
my $OddRt  = "[(\\bar s_\\alpha d_\\beta)_L(\\bar d_\\beta d_\\alpha)_R](x)";
my $OssR   = "[(\\bar s_\\alpha d_\\alpha)_L(\\bar s_\\beta s_\\beta)_R](x)";
my $OssRt  = "[(\\bar s_\\alpha d_\\beta)_L(\\bar u_\\beta u_\\alpha)_R](x)";
my $OccR   = "[(\\bar s_\\alpha d_\\alpha)_L(\\bar c_\\beta c_\\beta)_R](x)";
my $OccRt  = "[(\\bar s_\\alpha d_\\beta)_L(\\bar c_\\beta c_\\alpha)_R](x)";
my $OsdL   = "[\\bar s_\\alpha (1 - \\gamma_5) d_\\alpha](x)";
my $OsdR   = "[\\bar s_\\alpha (1 + \\gamma_5) d_\\alpha](x)";

# F(x) in tex format
#my @LFx =  ($OuuL,$OuuLt,$OddL,       $OssL,       $OccL,$OccLt,
#           $OuuR,$OuuRt,$OddR,$OddRt,$OssR,$OssRt,$OccR,$OccRt);
my @LFx =  ($OuuL,$OuuLt,$OddL,       $OssL,
            $OuuR,$OuuRt,$OddR,$OddRt,$OssR,$OssRt,
            $OsdL,$OsdR);

my $TuuL  = " ub_3^5 gmu2L^5^6 uv_3^6 db_4^7 gmu2L^7^8 sv_4^8";
my $TuuLt = " db_3^5 gmu2L^5^6 uv_3^6 ub_4^7 gmu2L^7^8 sv_4^8";
my $TddL  = " db_3^5 gmu2L^5^6 dv_3^6 db_4^7 gmu2L^7^8 sv_4^8";
my $TssL  = " sb_3^5 gmu2L^5^6 sv_3^6 db_4^7 gmu2L^7^8 sv_4^8";
my $TccL  = " cb_3^5 gmu2L^5^6 cv_3^6 db_4^7 gmu2L^7^8 sv_4^8";
my $TccLt = " db_3^5 gmu2L^5^6 cv_3^6 cb_4^7 gmu2L^7^8 sv_4^8";
my $TuuR  = " ub_3^5 gmu2R^5^6 uv_3^6 db_4^7 gmu2L^7^8 sv_4^8";
my $TuuRt = " - 2 db_3^5 gR^5^6 uv_3^6 ub_4^7 gL^7^8 sv_4^8";
my $TddR  = " db_3^5 gmu2R^5^6 dv_3^6 db_4^7 gmu2L^7^8 sv_4^8";
my $TddRt = " - 2 db_3^5 gR^5^6 dv_3^6 db_4^7 gL^7^8 sv_4^8";
my $TssR  = " sb_3^5 gmu2R^5^6 sv_3^6 db_4^7 gmu2L^7^8 sv_4^8";
my $TssRt = " - 2 db_3^5 gR^5^6 sv_3^6 sb_4^7 gL^7^8 sv_4^8";
my $TccR  = " cb_3^5 gmu2R^5^6 cv_3^6 db_4^7 gmu2L^7^8 sv_4^8";
my $TccRt = " - 2 db_3^5 gR^5^6 cv_3^6 cb_4^7 gL^7^8 sv_4^8";
my $TsdL  = " db_3^5 gR^5^6 sv_3^6";# ! L ↔ R?
my $TsdR  = " db_3^5 gL^5^6 sv_3^6";# ! L ↔ R?

# F(y)^\dag
my @Fy =  ($TuuL,$TuuLt,$TddL,       $TssL,       $TccL,$TccLt,
           $TuuR,$TuuRt,$TddR,$TddRt,$TssR,$TssRt,$TccR,$TccRt,
           $TsdL,$TsdR);
#my @Oy =  ($TuuL,$TuuLt,$TddL,       $TssL,
#           $TuuR,$TuuRt,$TddR,$TddRt,$TssR,$TssRt);
my @nameFy = ("uuLT","uuLmT","ddLT",        "ssLT",        "ccLT","ccLmT",
              "uuRT","uuRmT","ddRT","ddRmT","ssRT","ssRmT","ccRT","ccRmT",
              "blLT","blRT");

my $OTuuL  = "[(\\bar u_\\gamma u_\\gamma)_L(\\bar d_\\delta s_\\delta)_L](y)";
my $OTuuLt = "[(\\bar u_\\gamma u_\\delta)_L(\\bar d_\\delta s_\\gamma)_L](y)";
my $OTddL  = "[(\\bar d_\\gamma d_\\gamma)_L(\\bar d_\\delta s_\\delta)_L](y)";
my $OTssL  = "[(\\bar s_\\gamma s_\\gamma)_L(\\bar d_\\delta s_\\delta)_L](y)";
my $OTccL  = "[(\\bar c_\\gamma c_\\gamma)_L(\\bar d_\\delta s_\\delta)_L](y)";
my $OTccLt = "[(\\bar c_\\gamma c_\\delta)_L(\\bar d_\\delta s_\\gamma)_L](y)";
my $OTuuR  = "[(\\bar u_\\gamma u_\\gamma)_R(\\bar d_\\delta s_\\delta)_L](y)";
my $OTuuRt = "[(\\bar u_\\gamma u_\\delta)_R(\\bar d_\\delta s_\\gamma)_L](y)";
my $OTddR  = "[(\\bar d_\\gamma d_\\gamma)_R(\\bar d_\\delta s_\\delta)_L](y)";
my $OTddRt = "[(\\bar d_\\gamma d_\\delta)_R(\\bar d_\\delta s_\\gamma)_L](y)";
my $OTssR  = "[(\\bar s_\\gamma s_\\gamma)_R(\\bar d_\\delta s_\\delta)_L](y)";
my $OTssRt = "[(\\bar s_\\gamma s_\\delta)_R(\\bar d_\\delta s_\\gamma)_L](y)";
my $OTccR  = "[(\\bar c_\\gamma c_\\gamma)_R(\\bar d_\\delta s_\\delta)_L](y)";
my $OTccRt = "[(\\bar c_\\gamma c_\\delta)_R(\\bar d_\\delta s_\\gamma)_L](y)";
my $OTsdL  = "[\\bar s_\\gamma (1 + \\gamma_5) d_\\gamma](y)";# ! - ↔ +?
my $OTsdR  = "[\\bar s_\\gamma (1 - \\gamma_5) d_\\gamma](y)";# ! - ↔ +?

# F(y)^\dag in LaTeX format
my @LFy =  ($OTuuL,$OTuuLt,$OTddL,        $OTssL,        $OTccL,$OTccLt,
            $OTuuR,$OTuuRt,$OTddR,$OTddRt,$OTssR,$OTssRt,$OTccR,$OTccRt,
            $OTsdL,$OTsdR);
#my @Oy =  ($OTuuL,$OTuuLt,$OTddL,       $OTssL,
#           $OTuuR,$OTuuRt,$OTddR,$OTddRt,$OTssR,$OTssRt);


sub addifnew_acls {
  my $acls = shift;
  my @acls_objs = @_;

  my $itype = $acls->get_itype();
  my @tobjs = @{$acls_objs[$itype]};
  my $i = -1;
  #$acls->remove_conj_sign_flags();#
  foreach my $acls_i (@tobjs) {
    $i++;
    my $name_i = $acls_i->get_name_rep();
    my $name = $acls->get_name_rep();
    return @acls_objs if ( $name eq $name_i );
  }
  push (@tobjs,$acls);
  $acls_objs[$itype] = [@tobjs];
  return @acls_objs;
}

sub addifnew_tr {
  my $tr = shift;
  my @tr_objs = @_;
  my $itype = $tr->get_itype();
  my @trs = @{$tr_objs[$itype]};
  my $name = $tr->{name};
  for ( my $itr = 0 ; $itr <= $#trs ; $itr++ ) {
    my $name_i = $trs[$itr]->{name};
    return @tr_objs if ( $name eq $name_i );
  }
  push(@trs,$tr);
  #$trs[$#trs]->remove_conj_sign_flags();
  $tr_objs[$itype] = [@trs];
  return @tr_objs;
}

sub vector_formats {
  my @acls_objs = @_;
  foreach my $tobjs (@acls_objs) {
    foreach my $acls (@$tobjs) {
      foreach my $rep (@{$acls->{objs}}) {
        if ( $indiv_flag == 1 ) {
          for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
            print OUT "    std::vector<Rcomplex> " . "CON_" . $rep->{name}
            . $ixy . "(glb_v, Rcomplex(0.,0.));\n";
          }
        } else {
          print OUT "    std::vector<Rcomplex> " . "CON_" . $rep->{name}
                  . "(glb_v, Rcomplex(0.,0.));\n";
        }
      }
    }
  }
  return;
}

sub get_val_ind_munu {
  my $tr = shift;
  my $oddmu = $tr->get_n_oddmu("mu");
  my $oddnu = $tr->get_n_oddmu("nu");
  return "[munu]" if ( $oddmu == 1 && $oddnu == 1 );
  return "[mu]" if ( $oddmu == 1 && $oddnu == 0 );
  return "[nu]" if ( $oddmu == 0 && $oddnu == 1 );
  return "" if ( $oddmu == 0 && $oddnu == 0 );
  die;
}

sub get_xi_statement {
  my $rep = shift;
  my $ixy = shift;
  my $iyx = $ixy%2+1;

  # T4
  if ($rep->{name} =~ /c/) {
    foreach my $tr (@{$rep->{traces}}) {
      next unless ( $tr->{name} =~ /c/ );
      return " * conj(xi_l$ixy) * xi_l$iyx * conj(xi_c)"
          if ( $tr->{type} eq $tr_type_3y );
      last unless ( $tr->{type} eq $tr_type_1y || $tr->{type} eq $tr_type_4p );
      if ( $tr->{type} eq $tr_type_4p ) {
        return " * conj(xi_l$ixy) * xi_l$iyx * conj(xi_c)"
            if ( $tr->{flavors}->[2] eq "c" );
        # T2
        return " * xi_l$ixy * xi_l$iyx" if ( $tr->{flavors}->[0] eq "c" );
        die;
      }
      print "$tr->{name} $rep->{name}\n";
      die unless ( $tr->{type} eq $tr_type_1y );
      die unless ( $rep->{type} eq $rtype_31xy || $rep->{type} eq $rtype_211
                || $rep->{type} eq $rtype_21y );
      return " * conj(xi_l$ixy) * xi_l$iyx * conj(xi_c)" if ( $tr->{conjf} == 0 );
      return " * conj(xi_l$ixy) * xi_l$iyx * xi_c" if ( $tr->{conjf} == 1 );
    }
  }
  # T1
  return " * conj(xi_l$ixy) * xi_l$iyx"
      if ( $rep->{type} eq $rtype_2 || $rep->{type} eq $rtype_3x
        || $rep->{type} eq $rtype_21x );
  # T2
  return " * xi_l$ixy * xi_l$iyx"
      if ( $rep->{type} eq $rtype_4p || $rep->{type} eq $rtype_31yx
        || $rep->{type} eq $rtype_3y );
  if ( $rep->{type} eq $rtype_31xy || $rep->{type} eq $rtype_211
    || $rep->{type} eq $rtype_21y ) {
    foreach my $tr (@{$rep->{traces}}) {
      next unless ( $tr->{type} eq $tr_type_1y );
      return " * xi_l$ixy * xi_l$iyx" if ( $tr->{conjf} == 0 );
      return " * conj(xi_l$ixy) * conj(xi_l$iyx)" if ( $tr->{conjf} == 1 );
    }
  }
  # T3
  return " * xi_l$ixy * conj(xi_l$iyx)"
      if ( $rep->{type} eq $rtype_4 || $rep->{type} eq $rtype_22 );
  die;
}

sub add_products_munu {
  my @acls_objs = @_;
  foreach my $tobjs (@acls_objs) {
    foreach my $acls (@$tobjs) {
      foreach my $rep (@{$acls->{objs}}) {
        for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
          my $out = "";
          if ( $indiv_flag == 1 ) {
            $out = "              CON_" . $rep->{name} . $ixy . "[id] += ";
          } else {
            $out = "              CON_" . $rep->{name} . "[id] += ";
          }
          my @trs = @{$rep->{traces}};
#        my $tr1 = shift(@trs);
          my $flag = 0;
          foreach my $tr (@trs) {
            $out .= " * " unless ( $flag == 0 );
            $flag++;
            my $irand = $ixy;
            $irand = $irand%2+1 if ( $tr->{conjf} == 1 );
            $irand = "" if ( $tr->{name} =~ /F1x/ );
            $irand = "" if ( $tr->{name} =~ /F1y.+c$/ );
            if ( $tr->{conjf} == 1 ) {
              $out .= "conj($tr->{name}$irand" . &get_val_ind_munu($tr) . ")";
            } else {
              $out .= $tr->{name} . $irand . &get_val_ind_munu($tr);
            }
          }
          next unless ( $out =~ /munu/ );
          $out .= &get_xi_statement($rep,$ixy);
          print OUT "$out" . ";\n";
        }
      }
    }
  }
  return;
}

sub add_products_mu {
  my @acls_objs = @_;
  foreach my $tobjs (@acls_objs) {
    foreach my $acls (@$tobjs) {
      foreach my $rep (@{$acls->{objs}}) {
        for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
          my $out = "";
          if ( $indiv_flag == 1 ) {
            $out = "            CON_" . $rep->{name} . $ixy . "[id] += ";
          } else {
            $out = "            CON_" . $rep->{name} . "[id] += ";
          }
          my @trs = @{$rep->{traces}};
          my $flag = 0;
          foreach my $tr (@trs) {
            $out .= " * " unless ( $flag == 0 );
            $flag++;
            my $irand = $ixy;
            $irand = $irand%2+1 if ( $tr->{conjf} == 1 );
            $irand = "" if ( $tr->{name} =~ /F1x/ );
            $irand = "" if ( $tr->{name} =~ /F1y.+c$/ );
            if ( $tr->{conjf} == 1 ) {
              $out .= "conj($tr->{name}$irand" . &get_val_ind_munu($tr) . ")";
            } else {
              $out .= $tr->{name} . $irand . &get_val_ind_munu($tr);
            }
          }
          next if ( $out =~ /munu/ );
          next unless ( $out =~ /\[mu\]/ || $out =~ /\[nu\]/ );
          $out =~ s/\[nu\]/[mu]/g;
          $out .= &get_xi_statement($rep,$ixy);
          print OUT "$out" . ";\n";
        }
      }
    }
  }
  return;
}

sub add_products_scl {
  my @acls_objs = @_;
  foreach my $tobjs (@acls_objs) {
    foreach my $acls (@$tobjs) {
      foreach my $rep (@{$acls->{objs}}) {
        for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
          my $out = "";
          if ( $indiv_flag == 1 ) {
            $out = "          CON_" . $rep->{name} . $ixy . "[id] += ";
          } else {
            $out = "          CON_" . $rep->{name} . "[id] += ";
          }
          my @trs = @{$rep->{traces}};
          my $flag = 0;
          foreach my $tr (@trs) {
            $out .= " * " unless ( $flag == 0 );
            $flag++;
            my $irand = $ixy;
            $irand = $irand%2+1 if ( $tr->{conjf} == 1 );
            $irand = "" if ( $tr->{name} =~ /F1x/ );
            $irand = "" if ( $tr->{name} =~ /F1y.+c$/ );
            if ( $tr->{conjf} == 1 ) {
              $out .= "conj($tr->{name}$irand" . &get_val_ind_munu($tr) . ")";
            } else {
              $out .= $tr->{name} . $irand . &get_val_ind_munu($tr);
            }
          }
          next if ( $out =~ /munu/ );
          next if ( $out =~ /\[mu\]/ || $out =~ /\[nu\]/ );
          $out .= &get_xi_statement($rep,$ixy);
          print OUT "$out" . ";\n";
        }
      }
    }
  }
  return;
}

sub sum_over_all_nodes {
  my @acls_objs = @_;
  foreach my $tobjs (@acls_objs) {
    foreach my $acls (@$tobjs) {
      foreach my $rep (@{$acls->{objs}}) {
        if ( $indiv_flag == 1 ) {
          for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
            my $val = "CON_" . $rep->{name} . $ixy;
            print OUT "    glb_sum($val.data(), $val.size());\n";
          }
        } else {
          my $val = "CON_" . $rep->{name};
          print OUT "    glb_sum($val.data(), $val.size());\n";
        }
      }
    }
  }
  return;
}

sub output_on_binary_files {
  my @acls_objs = @_;
  my $icont = 0;
  foreach my $tobjs (@acls_objs) {
    foreach my $acls (@$tobjs) {
      foreach my $rep (@{$acls->{objs}}) {
        for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
          my $val = "CON_" . $rep->{name};
          $val .= $ixy if ( $indiv_flag == 1 );
          $val .= "[pt]";
          my $tmp = " $icont % nnode == nodeid ";
          print OUT "    if ($tmp) \{\n";
          print OUT "      char filename[1024];\n";
          print OUT "      sprintf(filename, \"Output/traj_%d_" . $rep->{name};
          print OUT $ixy if ( $indiv_flag == 1 );
          print OUT "\", meas_arg.TrajCur);\n";
          print OUT "      file.open(filename,std::ios::out | std::ios::binary);\n";
          print OUT "      for ( long long int pt = 0; pt < glb_v; pt++ ) \{\n";
          if ( $indiv_flag == 1 ) {
            print OUT "        $val /= nsrc;\n";
          } else {
            print OUT "        $val /= nsrc*2;\n";
          }
          print OUT "        auto temp = $val.real();\n";
          print OUT "        file.write((const char*)&temp, 8);\n";
          print OUT "      \}\n";
          print OUT "      file.close();\n";
          print OUT "      printf (\"Binary file saved on %s\\n\",filename);\n";
          print OUT "    \}\n";
          $icont++;
          last unless ( $indiv_flag == 1 );
        }
      }
    }
  }
  return;
}

sub traces_at_y {
  my @tr_objs = @_;
  my $itype = &Functions::inverse_array_char($tr_type_1y,[@tr_types]);
  my @trs = @{$tr_objs[$itype]};
  foreach my $tr (@trs) {
    my $G = $tr->{gammas}->[0];
    my $q = $tr->{flavors}->[0];
    my $nq = 2;
    $nq = 1 if ( $q =~ /c/ );
    for ( my $iq = 1 ; $iq <= $nq ; $iq++ ) {
      my $id = $iq;
      $id = "" if ( $q =~ /c/ );
      my $name = $tr->{name} . $id;
      if ( $G =~ /nu/ ) {
        print OUT "        std::vector<Rcomplex> " . $name
                . "(4, Rcomplex(0.,0.));\n";
        print OUT "        for ( nu = 0; nu <= 3; ++nu ) \{\n";
      } else {
        print OUT "        Rcomplex " . $name . "(0.,0.);\n";
        print OUT "        \{\n";
      }
      print OUT "          WilsonMatrix tmp = " . $q . "propyy" . $id . ";\n";
      if ( $G eq "gnuL" ) {
        print OUT "          tmp = tmp.glL(nu);// γν(1-γ5) S(y-y)\n";
        print OUT "          " . $name . "[nu] = Trace(tmp);\n";
      } elsif ( $G eq "gnuR" ) {
        print OUT "          tmp = tmp.glR(nu);// γν(1+γ5) S(y-y)\n";
        print OUT "          " . $name . "[nu] = Trace(tmp);\n";
      } elsif ( $G eq "gL" ) {
        print OUT "          tmp = tmp.glPL();// 0.5 (1-γ5) S(y-y)\n";
        print OUT "          tmp *= 2.0;\n";
        print OUT "          " . $name . " = Trace(tmp);\n";
      } elsif ( $G eq "gR" ) {
        print OUT "          tmp = tmp.glPR();// 0.5 (1+γ5) S(y-y)\n";
        print OUT "          tmp *= 2.0;\n";
        print OUT "          " . $name . " = Trace(tmp);\n";
      } else {
        print "$G\n";
        die;
      }
      print OUT "        \}\n";
    }
  }
  return;
}

sub declare_tr_at_x_andor_y {
  my $tr = shift;
  my $id = (shift||"");
  my $oddmu = $tr->get_n_oddmu("mu");
  my $oddnu = $tr->get_n_oddmu("nu");
  my $size = 4**($oddmu+$oddnu);
  if ( $size == 1 ) {
    print OUT "          Rcomplex " . $tr->{name} . $id . "(0.,0.);\n";
  } else {
    print OUT "          std::vector<Rcomplex> " . $tr->{name} . $id
            . "($size, Rcomplex(0.,0.));\n";
  }
  return ($oddmu,$oddnu);
}

sub print_mult_gamma_mu {
  my $G = shift;
  my $prop = shift;
  my $mu = shift;
  my $ind = shift;
  if ( $G eq ("g" . $mu . "L") ) {
    print OUT (" " x $ind) . "$prop = $prop.glL($mu);\n";
  } elsif ( $G eq ("g" . $mu . "R") ) {
    print OUT (" " x $ind) . "$prop = $prop.glR($mu);\n";
  } else {
    print "$G $mu\n";
    die;
  }
  return;
}

sub print_mult_gamma {
  my $G = shift;
  my $prop = shift;
  my $ind = shift;
  if ( $G eq "gL" ) {
    print OUT (" " x $ind) . "$prop = $prop.glPL();\n";
  } elsif ( $G eq "gR" ) {
    print OUT (" " x $ind) . "$prop = $prop.glPR();\n";
  }
  print OUT (" " x $ind) . "$prop *= 2.0;\n";
  return;
}

sub traces_at_x_and_y_4 {
  my @trs = @_;
  foreach my $tr (@trs) {
    # Tr[ propx1 propy1 propx2 propy2 ]
    # propx1 = Γ1 [lc]qpropxy[12]
    # propy1 = Γ2 [lc]qpropyx[21]
    # propx2 = Γ3 [lc]qpropxy[12]
    # propy2 = Γ4 [lc]qpropyx[21]
    my @G = @{$tr->{gammas}};
    my @F = @{$tr->{flavors}};
    for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
      my $iyx = $ixy%2+1;
      print OUT "          Rcomplex " . $tr->{name} . $ixy . "(0.,0.);\n";
      unless ( $G[0] =~ /mu/ ) {
        print OUT "          propx1 = " . $F[0] . "propxy" . $ixy . ";\n";
        &print_mult_gamma($G[0],"propx1",10);
        print OUT "          propx2 = " . $F[2] . "propxy" . $ixy . ";\n";
        &print_mult_gamma($G[2],"propx2",10);
      }
      unless ( $G[1] =~ /nu/ ) {
        print OUT "          propy1 = " . $F[1] . "propyx" . $iyx . ";\n";
        &print_mult_gamma($G[1],"propy1",10);
        print OUT "          propy2 = " . $F[3] . "propyx" . $iyx . ";\n";
        &print_mult_gamma($G[3],"propy2",10);
      }
      if ( $G[0] =~ /mu/ ) {
        print OUT "          for ( mu = 0; mu <= 3; ++mu ) \{\n";
        print OUT "            propx1 = " . $F[0] . "propxy" . $ixy . ";\n";
        &print_mult_gamma_mu($G[0],"propx1","mu",12);
        print OUT "            propx2 = " . $F[2] . "propxy" . $ixy . ";\n";
        &print_mult_gamma_mu($G[2],"propx2","mu",12);
        if ( $G[1] =~ /nu/ ) {
          print OUT "            for ( nu = 0; nu <= 3; ++nu ) \{\n";
          print OUT "              propy1 = " . $F[1] . "propyx" . $iyx . ";\n";
          &print_mult_gamma_mu($G[1],"propy1","nu",14);
          print OUT "              propy2 = " . $F[3] . "propyx" . $iyx . ";\n";
          &print_mult_gamma_mu($G[3],"propy2","nu",14);
          print OUT "              " . $tr->{name} . $ixy
                  . " += get_F4(propx1,propy1,propx2,propy2);\n";
          print OUT "            \}// nu\n";
        } else {
          print OUT "            " . $tr->{name} . $ixy
                  . " += get_F4(propx1,propy1,propx2,propy2);\n";
        }
        print OUT "          \}// mu\n";
      } elsif ( $G[1] =~ /nu/ ) {
        print OUT "          for ( nu = 0; nu <= 3; ++nu ) \{\n";
        print OUT "            propy1 = " . $F[1] . "propyx" . $iyx . ";\n";
        &print_mult_gamma_mu($G[1],"propy1","nu",12);
        print OUT "            propy2 = " . $F[3] . "propyx" . $iyx . ";\n";
        &print_mult_gamma_mu($G[3],"propy2","nu",12);
        print OUT "            " . $tr->{name} . $ixy
                . " += get_F4(propx1,propy1,propx2,propy2);\n";
        print OUT "          \}// nu\n";
      } else {
        print OUT "          " . $tr->{name} . $ixy
                . " = get_F4(propx1,propy1,propx2,propy2);\n";
      }
    }
  }
  return;
}

sub traces_at_x_and_y_4p {
  my @trs = @_;
  foreach my $tr (@trs) {
    my @G = @{$tr->{gammas}};
    my @F = @{$tr->{flavors}};
    # Tr[ propx1 propx2 propy1 propy2 ]
    # propx1 = Γ1 [lc]qpropxx
    # propx2 = Γ2 [lc]qpropxy
    # propy1 = Γ3 [lc]qpropyy
    # propy2 = Γ4 [lc]qpropyx
    for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
      my $iyx = $ixy%2+1;
      my $iyy = $ixy;
      $iyy = "" if ( $F[2] =~ /c/ );
      print OUT "          Rcomplex " . $tr->{name} . $ixy . "(0.,0.);\n";
      unless ( $G[0] =~ /mu/ ) {
        print OUT "          propx1 = " . $F[0] . "propxx;\n";
        &print_mult_gamma($G[0],"propx1",10);
        print OUT "          propx2 = " . $F[1] . "propxy" . $ixy . ";\n";
        &print_mult_gamma($G[1],"propx2",10);
      }
      unless ( $G[2] =~ /nu/ ) {
        print OUT "          propy1 = " . $F[2] . "propyy" . $iyy . ";\n";
        &print_mult_gamma($G[2],"propy1",10);
        print OUT "          propy2 = " . $F[3] . "propyx" . $iyx . ";\n";
        &print_mult_gamma($G[3],"propy2",10);
      }
      if ( $G[0] =~ /mu/ ) {
        print OUT "          for ( mu = 0; mu <= 3; ++mu ) \{\n";
        print OUT "            propx1 = " . $F[0] . "propxx;\n";
        &print_mult_gamma_mu($G[0],"propx1","mu",12);
        print OUT "            propx2 = " . $F[1] . "propxy" . $ixy . ";\n";
        &print_mult_gamma_mu($G[1],"propx2","mu",12);
        if ( $G[2] =~ /nu/ ) {
          print OUT "            for ( nu = 0; nu <= 3; ++nu ) \{\n";
          print OUT "              propy1 = " . $F[2] . "propyy" . $iyy . ";\n";
          &print_mult_gamma_mu($G[2],"propy1","nu",14);
          print OUT "              propy2 = " . $F[3] . "propyx" . $iyx . ";\n";
          &print_mult_gamma_mu($G[3],"propy2","nu",14);
          print OUT "              " . $tr->{name} . $ixy
                  . " += get_F4(propx1,propx2,propy1,propy2);\n";
          print OUT "            \}// nu\n";
        } else {
          print OUT "            " . $tr->{name} . $ixy
                  . " += get_F4(propx1,propx2,propy1,propy2);\n";
        }
        print OUT "          \}// mu\n";
      } elsif ( $G[2] =~ /nu/ ) {
        print OUT "          for ( nu = 0; nu <= 3; ++nu ) \{\n";
        print OUT "            propy1 = " . $F[2] . "propyy" . $iyy . ";\n";
        &print_mult_gamma_mu($G[2],"propy1","nu",12);
        print OUT "            propy2 = " . $F[3] . "propyx" . $iyx . ";\n";
        &print_mult_gamma_mu($G[3],"propy2","nu",12);
        print OUT "            " . $tr->{name} . $ixy
                . " += get_F4(propx1,propx2,propy1,propy2);\n";
        print OUT "          \}// nu\n";
      } else {
        print OUT "          " . $tr->{name} . $ixy
                . " = get_F4(propx1,propx2,propy1,propy2);\n";
      }
    }
  }
  return;
}

sub traces_at_x_and_y_3x {
  my @trs = @_;
  foreach my $tr (@trs) {
    my @G = @{$tr->{gammas}};
    my @F = @{$tr->{flavors}};
    # Tr[ propx1 propx2 propy1 ]
    # propx1 = Γ1 [lc]qpropxx
    # propx2 = Γ2 [lc]qpropxy
    # propy1 = Γ3 [lc]qpropyx
    for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
      my $iyx = $ixy%2+1;
      my @oddmunu = &declare_tr_at_x_andor_y($tr,$ixy);
      unless ( $G[0] =~ /mu/ ) {
        print OUT "          propx1 = " . $F[0] . "propxx;\n";
        &print_mult_gamma($G[0],"propx1",10);
        print OUT "          propx2 = " . $F[1] . "propxy" . $ixy . ";\n";
        &print_mult_gamma($G[1],"propx2",10);
      }
      unless ( $G[2] =~ /nu/ ) {
        print OUT "          propy1 = " . $F[2] . "propyx" . $iyx . ";\n";
        &print_mult_gamma($G[2],"propy1",10);
      }
      if ( $G[0] =~ /mu/ ) {
        print OUT "          for ( mu = 0; mu <= 3; ++mu ) \{\n";
        print OUT "            propx1 = " . $F[0] . "propxx;\n";
        &print_mult_gamma_mu($G[0],"propx1","mu",12);
        print OUT "            propx2 = " . $F[1] . "propxy" . $ixy . ";\n";
        &print_mult_gamma_mu($G[1],"propx2","mu",12);
        if ( $G[2] =~ /nu/ ) {
          print OUT "            for ( nu = 0; nu <= 3; ++nu ) \{\n";
          print OUT "              propy1 = " . $F[2] . "propyx" . $iyx . ";\n";
          &print_mult_gamma_mu($G[2],"propy1","nu",14);
          print OUT "              " . $tr->{name} . $ixy
                  . "[nu] += get_F3(propx1,propx2,propy1);\n";
          print OUT "            \}// nu\n";
        } else {
          print OUT "            " . $tr->{name} . $ixy
                  . " += get_F3(propx1,propx2,propy1);\n";
        }
        print OUT "          \}// mu\n";
      } elsif ( $G[2] =~ /nu/ ) {
        print OUT "          for ( nu = 0; nu <= 3; ++nu ) \{\n";
        print OUT "            propy1 = " . $F[2] . "propyx" . $iyx . ";\n";
        &print_mult_gamma_mu($G[2],"propy1","nu",12);
        print OUT "            " . $tr->{name} . $ixy
                . "[nu] += get_F3(propx1,propx2,propy1);\n";
        print OUT "          \}// nu\n";
      } else {
        print OUT "          " . $tr->{name} . $ixy
                . " = get_F3(propx1,propx2,propy1);\n";
      }
    }
  }
  return;
}

sub traces_at_x_and_y_3y {
  my @trs = @_;
  foreach my $tr (@trs) {
    my @G = @{$tr->{gammas}};
    my @F = @{$tr->{flavors}};
    # Tr[ propy1 propy2 propx1 ]
    # propy1 = Γ1 [lc]qpropyy
    # propy2 = Γ2 [lc]qpropyx
    # propx1 = Γ3 [lc]qpropxy
    for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
      my $iyx = $ixy%2+1;
      my $iyy = $ixy;
      $iyy = "" if ( $F[0] =~ /c/ );
      my @oddmunu = &declare_tr_at_x_andor_y($tr,$ixy);
      unless ( $G[0] =~ /nu/ ) {
        print OUT "          propy1 = " . $F[0] . "propyy" . $iyy . ";\n";
        &print_mult_gamma($G[0],"propy1",10);
        print OUT "          propy2 = " . $F[1] . "propyx" . $iyx . ";\n";
        &print_mult_gamma($G[1],"propy2",10);
      }
      unless ( $G[2] =~ /mu/ ) {
        print OUT "          propx1 = " . $F[2] . "propxy" . $ixy . ";\n";
        &print_mult_gamma($G[2],"propx1",10);
      }
      if ( $G[0] =~ /nu/ ) {
        print OUT "          for ( nu = 0; nu <= 3; ++nu ) \{\n";
        print OUT "            propy1 = " . $F[0] . "propyy" . $iyy . ";\n";
        &print_mult_gamma_mu($G[0],"propy1","nu",12);
        print OUT "            propy2 = " . $F[1] . "propyx" . $iyx . ";\n";
        &print_mult_gamma_mu($G[1],"propy2","nu",12);
        if ( $G[2] =~ /mu/ ) {
          print OUT "            for ( mu = 0; mu <= 3; ++mu ) \{\n";
          print OUT "              propx1 = " . $F[2] . "propxy" . $ixy . ";\n";
          &print_mult_gamma_mu($G[2],"propx1","mu",14);
          print OUT "              " . $tr->{name} . $ixy
                  . "[mu] += get_F3(propy1,propy2,propx1);\n";
          print OUT "            \}// mu\n";
        } else {
          print OUT "            " . $tr->{name} . $ixy
                  . " += get_F3(propy1,propy2,propx1);\n";
        }
        print OUT "          \}// nu\n";
      } elsif ( $G[2] =~ /mu/ ) {
        print OUT "          for ( mu = 0; mu <= 3; ++mu ) \{\n";
        print OUT "            propx1 = " . $F[2] . "propxy" . $ixy . ";\n";
        &print_mult_gamma_mu($G[2],"propx1","mu",12);
        print OUT "            " . $tr->{name} . $ixy
                . "[mu] += get_F3(propy1,propy2,propx1);\n";
        print OUT "          \}// mu\n";
      } else {
        print OUT "          " . $tr->{name} . $ixy
                . " = get_F3(propy1,propy2,propx1);\n";
      }
    }
  }
  return;
}

sub traces_at_x_and_y_2 {
  my @trs = @_;
  foreach my $tr (@trs) {
    my @G = @{$tr->{gammas}};
    my @F = @{$tr->{flavors}};
    # Tr[ propx1 propy1 ]
    # propx1 = Γ1 [lc]qpropxy
    # propy1 = Γ2 [lc]qpropyx
    for ( my $ixy = 1 ; $ixy <= 2 ; $ixy++ ) {
      my $iyx = $ixy%2+1;
      my @oddmunu = &declare_tr_at_x_andor_y($tr,$ixy);
      unless ( $G[0] =~ /mu/ ) {
        print OUT "          propx1 = " . $F[0] . "propxy" . $ixy . ";\n";
        &print_mult_gamma($G[0],"propx1",10);
      }
      unless ( $G[1] =~ /nu/ ) {
        print OUT "          propy1 = " . $F[1] . "propyx" . $iyx . ";\n";
        &print_mult_gamma($G[1],"propy1",10);
      }
      if ( $G[0] =~ /mu/ ) {
        print OUT "          for ( mu = 0; mu <= 3; ++mu ) \{\n";
        print OUT "            propx1 = " . $F[0] . "propxy" . $ixy . ";\n";
        &print_mult_gamma_mu($G[0],"propx1","mu",12);
        if ( $G[1] =~ /nu/ ) {
          print OUT "            for ( nu = 0; nu <= 3; ++nu ) \{\n";
          print OUT "              propy1 = " . $F[1] . "propyx" . $iyx . ";\n";
          &print_mult_gamma_mu($G[1],"propy1","nu",14);
          print OUT "              int munu = 4*mu + nu;\n";
          print OUT "              " . $tr->{name} . $ixy
                  . "[munu] += Trace(propx1,propy1);\n";
          print OUT "            \}// nu\n";
        } else {
          print OUT "            " . $tr->{name} . $ixy
                  . "[mu] += Trace(propx1,propy1);\n";
        }
        print OUT "          \}// mu\n";
      } elsif ( $G[1] =~ /nu/ ) {
        print OUT "          for ( nu = 0; nu <= 3; ++nu ) \{\n";
        print OUT "            propy1 = " . $F[1] . "propyx" . $iyx . ";\n";
        &print_mult_gamma_mu($G[1],"propy1","nu",12);
        print OUT "            " . $tr->{name} . $ixy
                . "[nu] += Trace(propx1,propy1);\n";
        print OUT "          \}// nu\n";
      } else {
        print OUT "          " . $tr->{name} . $ixy
                . " = Trace(propx1,propy1);\n";
      }
    }
  }
  return;
}

sub traces_at_x_1x {
  my @trs = @_;
  foreach my $tr (@trs) {
    my @oddmunu = &declare_tr_at_x_andor_y($tr);
    # Tr[ propx1 ]
    # propx1 = Γ1 [lc]qpropxx
    my @G = @{$tr->{gammas}};
    my @F = @{$tr->{flavors}};
    unless ( $G[0] =~ /mu/ ) {
      print OUT "          propx1 = " . $F[0] . "propxx;\n";
      &print_mult_gamma($G[0],"propx1",10);
    }
    if ( $G[0] =~ /mu/ ) {
      print OUT "          for ( mu = 0; mu <= 3; ++mu ) \{\n";
      print OUT "            propx1 = " . $F[0] . "propxx;\n";
      &print_mult_gamma_mu($G[0],"propx1","mu",12);
      print OUT "            " . $tr->{name}
              . "[mu] += Trace(propx1);\n";
      print OUT "          \}// mu\n";
    } else {
      print OUT "          " . $tr->{name} . " = Trace(propx1);\n";
    }
  }
  return;
}

sub traces_at_x_andor_y {
  my @tr_objs = @_;
  my $itype_4  = &Functions::inverse_array_char($tr_type_4, [@tr_types]);
  my $itype_4p = &Functions::inverse_array_char($tr_type_4p,[@tr_types]);
  my $itype_3x = &Functions::inverse_array_char($tr_type_3x,[@tr_types]);
  my $itype_3y = &Functions::inverse_array_char($tr_type_3y,[@tr_types]);
  my $itype_2  = &Functions::inverse_array_char($tr_type_2, [@tr_types]);
  my $itype_1x = &Functions::inverse_array_char($tr_type_1x,[@tr_types]);
  my $itype_1y = &Functions::inverse_array_char($tr_type_1y,[@tr_types]);
  for ( my $itype = 0 ; $itype <= $#tr_objs ; $itype++ ) {
    next if ( $itype == $itype_1y );
    my @trs = @{$tr_objs[$itype]};
    &traces_at_x_and_y_4(@trs)  if ( $itype == $itype_4 );
    &traces_at_x_and_y_4p(@trs) if ( $itype == $itype_4p );
    &traces_at_x_and_y_3x(@trs) if ( $itype == $itype_3x );
    &traces_at_x_and_y_3y(@trs) if ( $itype == $itype_3y );
    &traces_at_x_and_y_2(@trs)  if ( $itype == $itype_2 );
    &traces_at_x_1x(@trs) if ( $itype == $itype_1x );
  }
  return;
}

my @corr_objs = ();# will be used to create analysis code
my @acls_objs = &Functions::new_ref_array(7);# will be used for tex and measurement codes
my $out_fname = ( shift || "contractions_auto.tex" );
open(OUT,"> $out_fname") or die "Failed to open $out_fname\n";
print OUT "\\begin{align}\n";
for ( my $c = 0 ; $c <= $#Fx ; $c++ ) {
  my $O1 = $Fx[$c];
  my $nameO1 = $nameFx[$c];
  my $flag = 1;
  my @objs_c = ();
  for ( my $l = 0 ; $l <= $#Fy ; $l++ ) {
    my $O2 = $Fy[$l];
    my $nameO2 = $nameFy[$l];
    print "########################################################\n";
    print ( ($c+1) . "-" . ($l+1) . "\n" );
    my $Opr1 = $O1;
    my $Opr2 = $O2;
    $Opr1 =~ s/_\d//g;
    $Opr1 =~ s/\^\d//g;
    $Opr2 =~ s/_\d//g;
    $Opr2 =~ s/\^\d//g;
    print "O1: $Opr1\n";
    print "O2: $Opr2\n";

    my $corr = CorrManage->new($O1,$O2);
    $corr->show_object();

    my @cont = $corr->trace_contraction(2);
    next if ( $#cont == -1 );

    if ( $c == 0 && $l == 0 ) {
    } elsif ( $flag == 1 ) {
      $flag = 0;
      print OUT "\\\\[3mm]\n";
    } else {
      print OUT "\\\\\n";
#      print OUT "\\\\*\n";
    }

#    print OUT ( "%" . ($c+1) . "-" . ($l+1) . "\n" );
    print OUT ( "\\mbox\{(" . ($c+1) . "-" . ($l+1) . "):\}\\ \n" );
    print OUT "\\Big\\langle\n";
    print OUT "$LFx[$c]&\n";
    print OUT "$LFy[$l]\\big|_{(I)}\n";
    print OUT "\\Big\\rangle\n";
    print OUT "\\notag\\\\*\n";
    print OUT "&=";

    my @aclss = ();
    foreach my $con(@cont) {
      my $coef = 0;
      $coef = $1 if ( $con =~ /^([^T]+)T/ );
      my $acls = FormatManage->new($con);
      my $itype = $acls->get_itype();
      push(@aclss,[$itype,$acls,$coef]);
      $acls->get_sign_rep();
      my $acls2 = FormatManage->copy_general($acls);
      @acls_objs = &addifnew_acls($acls2,@acls_objs);
      $acls->get_sign_rep();
    }# foreach $con

    # Checked
    @aclss = sort {$a->[0] <=> $b->[0]} @aclss;
    my @texform = ();

    my $decomp_fname = $decomp_dir . $nameO1 . "_" . $nameO2;
    open(DCP,"> $decomp_fname") or die "Failed to open DCP: $decomp_fname\n";
#    foreach my $tex (@texf) {
    for ( my $it = 0 ; $it <= $#aclss ; $it++ ) {
      my $acls = $aclss[$it]->[1];
      $acls->get_sign_rep();
      my $tex = $acls->get_texf_rep();
      my $coef = $aclss[$it]->[2];
      $coef = 1 if ( $coef =~ /^\s*\+\s*$/ );
      $coef = -1 if ( $coef =~ /^\s*\-\s*$/ );
      $coef *= $acls->get_sign_rep();
      $acls->{coef} = $coef;
      print DCP "$coef " . $acls->get_name_rep() . "\n";
      if ( $coef == 1 ) {
        $coef = " + " ;
      } elsif ( $coef == -1 ) {
        $coef = " - ";
      } elsif ( $coef > 0 ) {
        $coef = " + $coef";
      }
      $tex = $coef . $tex;
      print OUT "\n";
      print OUT "\\notag\\\\*&\\ \\ \n" unless ( $it == 0 );
      print OUT "$tex";
      $aclss[$it] = $acls;
    }
    close(DCP);
    print "Decomposition written on $decomp_fname\n";
#    $corr->{texf} = [@aclss];
    print OUT ",\n" unless ( $c == $#Fx && $l == $#Fy );
    print OUT ".\n" if ( $c == $#Fx && $l == $#Fy );
    push (@objs_c,[@aclss]);
#    exit 0;
  }
  push @corr_objs, [@objs_c];
#  exit 0;
}
print OUT "\\end{align}\n\n";

print "\n";
print "=  Outputting needed contractions\n";
print "&  Extracting of needed traces\n";
my @tr_objs = &Functions::new_ref_array(6);# will be used for tex and measurement codes
for ( my $itype = 0 ; $itype <= $#acls_types ; $itype++ ) {
  my @tobjs = @{$acls_objs[$itype]};
  my $nindiv = $#tobjs+1;
  my $ntot = 0;
  my $nc = 0;
  foreach my $acls (@tobjs) {
    $ntot += $acls->{nobjs};
    $nc++ unless ( $acls->{nc} == 0 );
  }
  print OUT $acls_types[$itype] . "-type: $nindiv independent contractions (including $nc charm-contained ones) / $ntot total contractions:\n";
  print OUT "\\begin{align*}\n";
  foreach my $acls (@tobjs) {
    my @reps = @{$acls->{objs}};
    if ( $#reps == 0 ) {
      print OUT "&" . $reps[0]->{texf} . "\\\\\n";
    } elsif ( $#reps == 1 ) {
      print OUT "&" . $reps[0]->{texf} . "\n";
      print OUT "=" . $reps[1]->{texf} . "\\\\\n";
      print "$reps[0]->{name} \\= $reps[1]->{name}\n";
    } else {
      print "$#reps\n";
      die;
    }
    foreach my $rep (@reps) {
      foreach my $tr (@{$rep->{traces}}) {
        @tr_objs = &addifnew_tr($tr,@tr_objs);
      }
    }
  }
  print OUT "\\end{align*}\n";
}
print OUT "\n";

for ( my $itype = 0 ; $itype <= $#tr_ttypes ; $itype++ ) {
  my @trs = @{$tr_objs[$itype]};
  my $nindiv = $#trs + 1;
  print OUT $tr_ttypes[$itype] . "-type: $nindiv traces to be taken:\n";
  print OUT "\\begin{align*}\n";
  foreach my $tr (@trs) {
    print OUT "&" . $tr->{texf} . ",\\\\\n";
  }
  print OUT "\\end{align*}\n";
}
close (OUT);
print "\n";
print "Linear combination of flavor-basis correlators saved in $decomp_dir\n";
print "Writing on TeX file $out_fname completed!\n";

open(IN,"< $tmp_meas_code") or die "Failed to open $tmp_meas_code\n";
open(OUT,"> $meas_code") or die "Failed to open $meas_code\n";
while(<IN>) {
  s/__NR__/$nr/;
  s/__NSRC_S__/$nsrc_s/;
  s/__NSRC_T__/$nsrc_t/;
  print OUT $_;
  &vector_formats(@acls_objs) if ( /__VECTOR_FORMATS__/ );
  &traces_at_y(@tr_objs) if ( /__TRACES_AT_Y__/ );
  &traces_at_x_andor_y(@tr_objs) if ( /__TRACES_AT_X_ANDOR_Y__/ );
  &add_products_munu(@acls_objs)
      if ( /__ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_MUNU__/ );
  &add_products_mu(@acls_objs)
      if ( /__ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_MU__/ );
  &add_products_scl(@acls_objs)
      if ( /__ADD_PRODUCTS_OF_TRACES_TO_CONTRACTIONS_SCL__/ );
  &sum_over_all_nodes(@acls_objs) if ( /__SUM_OVER_ALL_NODES__/ );
  &output_on_binary_files(@acls_objs) if ( /__OUTPUT_ON_BINARY_FILES__/ );
}
close(OUT);
close(IN);
print "Measurement code $meas_code created!\n";

open (OUT,"> $ch_list") or die "Failed to open $ch_list\n";
foreach my $tobjs (@acls_objs) {
  foreach my $acls (@$tobjs) {
    foreach my $rep (@{$acls->{objs}}) {
      # if ( $indiv_flag == 1 ) {
      # for (... ) ...
      # next;
      # }
      print OUT " " . $rep->{name};
    }
    print OUT "\n";
  }
}
close(OUT);
print "List of Channels written on $ch_list\n";

exit 0;
