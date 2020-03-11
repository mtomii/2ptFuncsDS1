package TraceManage;
use strict;
use warnings;

our $write_flag = 4;
sub get_write_flag { return $write_flag };

our $type_4   = "F4";
our $type_4p  = "F4p";
our $type_3x  = "F3x";
our $type_3y  = "F3y";
our $type_2   = "F2";
our $type_1x  = "F1x";
our $type_1y  = "F1y";
sub get_type_4 { return $type_4 }
sub get_type_4p { return $type_4p }
sub get_type_3x { return $type_3x }
sub get_type_3y { return $type_3y }
sub get_type_2 { return $type_2 }
sub get_type_1x { return $type_1x }
sub get_type_1y { return $type_1y }
our @types = ($type_4,$type_4p,$type_3x,$type_3y,
              $type_2,$type_1x,$type_1y);
sub get_types { return @types }

our $ttype_4   = "\$F_4\$";
our $ttype_4p  = "\$F_4'\$";
our $ttype_3x  = "\$F_{3,x}\$";
our $ttype_3y  = "\$F_{3,y}\$";
our $ttype_2   = "\$F_2\$";
our $ttype_1x  = "\$F_{1,x}\$";
our $ttype_1y  = "\$F_{1,y}\$";
#sub get_ttype_4 { return $ttype_4 }
#sub get_ttype_4p { return $ttype_4p }
#sub get_ttype_3x { return $ttype_3x }
#sub get_ttype_3y { return $ttype_3y }
#sub get_ttype_2 { return $ttype_2 }
#sub get_ttype_1x { return $ttype_1x }
sub get_ttype_1y { return $ttype_1y }
our @ttypes = ($ttype_4,$ttype_4p,$ttype_3x,$ttype_3y,
               $ttype_2,$ttype_1x,$ttype_1y);
sub get_ttypes { return @ttypes }

sub new {
  my $class = shift;
  my $type = shift;
  my @gamfl = @_;
  if ( $write_flag > 1 ) {
    print "### TraceManage->new with rep type $type\n";
    my @G = @{$gamfl[0]};
    print "     gammas:  @G\n";
    my @F = @{$gamfl[1]};
    print "     flavors: @F\n";
  }
  my $this = {
               type => $type,
               gammas => $gamfl[0],# may be changed by $this->set_rep_gamfl();
               flavors => $gamfl[1],# may be changed by $this->set_rep_gamfl();
               texf => undef,
               name => undef,
               conjf => 0,# may be changed by $this->set_rep_gamfl();
               sign => 1
  };

  bless $this, $class;

  $this->set_rep_gamfl();
  $this->set_texf();
  $this->set_name();
  if ( $write_flag > 1 ) {
    print "### texf of Trace: $this->{texf}\n";
    print "### name of Trace: $this->{name}\n";
  }

  return $this;
}

sub get_iq {
  my $F = shift;
  return 1 if ( $F eq "l" );
  return 2 if ( $F eq "s" );
  return 3 if ( $F eq "c" );
  print "$F\n";
  die;
}

sub get_ig {
  my $G = shift;
  return 3 if ( $G eq "gL" );
  return 4 if ( $G eq "gR" );
  return 1 if ( $G =~ /^gmu\d*L$/ || $G =~ /^gnuL$/ );
  return 2 if ( $G =~ /^gmu\d*R$/ || $G =~ /^gnuR$/ );
  print "$G\n";
  die;
}

sub choose_rep_gamfl {
  my $this = shift;
  my $Gref = shift;
  my $Fref = shift;
  my $conjf = shift;

  my @F  = @$Fref;
  my @F0 = @{$this->{flavors}};
  for ( my $i = 0 ; $i <= $#F ; $i++ ) {
    next if ( $F[$i] eq $F0[$i] );
    my $iq  = &get_iq($F[$i]);
    my $iq0 = &get_iq($F0[$i]);
    return if ( $iq0 < $iq );
    $this->{flavors} = $Fref;
    $this->{gammas}  = $Gref;
    if ( $write_flag > 2 ) {
      print "### Update:\n";
      my @Fnew = @$Fref;
      my @Gnew = @$Gref;
      print "     gammas:  @Gnew\n";
      print "     flavors: @Fnew\n";
    }
    #return if ( $this->{conjf} == -1 );
    $this->{conjf} = $conjf;
    print "     conjugate flag: $conjf\n" if ( $write_flag > 2 );
    if ( $conjf == 1 ) {
      $this->{sign} = 1;
      foreach my $gam (@$Gref) {
        $this->{sign} *= -1 if ( $gam =~ /mu/ );
      }
      print "     Overall sign: $this->{sign}\n" if ( $write_flag > 2 );
    }
    return;
  }

  my @G  = @$Gref;
  my @G0 = @{$this->{gammas}};
  for ( my $i = 0 ; $i <= $#G ; $i++ ) {
    next if ( $G[$i] eq $G0[$i] );
    my $ig  = &get_ig($G[$i]);
    my $ig0 = &get_ig($G0[$i]);
    return if ( $ig0 < $ig );
    $this->{flavors} = $Fref;
    $this->{gammas}  = $Gref;
    if ( $write_flag > 2 ) {
      print "### Update:\n";
      my @Fnew = @$Fref;
      my @Gnew = @$Gref;
      print "     gammas:  @Gnew\n";
      print "     flavors: @Fnew\n";
    }
    #return if ( $this->{conjf} == -1 );
    $this->{conjf} = $conjf;
    print "     conjugate flag: $conjf\n" if ( $write_flag > 2 );
    if ( $conjf == 1 ) {
      $this->{sign} = 1;
      #my @test = @$Gref;
      #print "######  ";
      foreach my $gam (@$Gref) {
        #print "  $gam";
        $this->{sign} *= -1 if ( $gam =~ /mu/ );
      }
      #print "\n";
      print "     Overall sign: $this->{sign}\n" if ( $write_flag > 2 );
    }
    return;
  }
  #$this->{conjf} = -1 if ( $this->{conjf} == 0 && $conjf == 1 ); to indicate it is a real value
  return;
}

sub get_g5Hconjs {
  my @G = @_;
  my @result = ();
  foreach my $Gam(@G) {
    # g[LR]: γ5 [(1 ± γ5)]† γ5 = γμ (1 ± γ5)
    # gmu[12][LR]: γ5 [γμ (1 ± γ5)]† γ5 = - γμ (1 ∓ γ5)
    # The sign '-' is taken into account in sub choose_rep_gamfl
    my $Gam1 = $Gam;
    if ( $Gam =~ /mu/ ) {
      if ( $Gam =~ /L$/ ) {
        $Gam1 =~ s/L$/R/;
      } elsif ( $Gam =~ /R$/ ) {
        $Gam1 =~ s/R$/L/;
      } else {
        print "$Gam\n";
        die;
      }
    } else {
      die unless ( $Gam eq "gL" || $Gam eq "gR" );
    }
    push(@result,$Gam1);
  }
  return @result;
}

sub tex_gamma {
  my $gamma = shift;
  return "\\Gamma_\\mu^-" if ( $gamma eq "gmuL" );
  return "\\Gamma_\\mu^+" if ( $gamma eq "gmuR" );
  return "\\Gamma_\\nu^-" if ( $gamma eq "gnuL" );
  return "\\Gamma_\\nu^+" if ( $gamma eq "gnuR" );
  return "1-\\gamma_5" if ( $gamma eq "gL" );
  return "1+\\gamma_5" if ( $gamma eq "gR" );
  print "$gamma\n";
  die;
}

##     ##
#  F_4  #
##     ##

sub get_cand_star_4 {
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) Γ4 Sf4(y-x) ]
  # This sub translates:
  #    (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4) -> (Γ1★,Γ4★,Γ3★,Γ2★)(F4,F3,F2,F1)
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  @G = &get_g5Hconjs(@G);
  @G = ($G[0],$G[3],$G[2],$G[1]);
  @F = ($F[3],$F[2],$F[1],$F[0]);
  return([@G],[@F]);
}

sub clean_munu_4 {
  my $this = shift;
  my @G = @{$this->{gammas}};
  $G[0] =~ s/mu\d/mu/;
  $G[1] =~ s/mu\d/nu/;
  $G[2] =~ s/mu\d/mu/;
  $G[3] =~ s/mu\d/nu/;
  $this->{gammas} = [@G];
  return;
}

sub set_rep_gamfl_4 {
  my $this = shift;
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) Γ4 Sf4(y-x) ]
  #   Cand:  (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4), (Γ3,Γ4,Γ1,Γ2)(F3,F4,F1,F2)
  #   Cand*: (Γ1★,Γ4★,Γ3★,Γ2★)(F4,F3,F2,F1), ...
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  my @cand_1 = ([@G],[@F]);
  my @cand_2 = ( [$G[2],$G[3],$G[0],$G[1]], [$F[2],$F[3],$F[0],$F[1]] );
  my @cand_3 = &get_cand_star_4(@cand_1);
  my @cand_4 = &get_cand_star_4(@cand_2);

  # Update to representative gamfl with conjf
  $this->choose_rep_gamfl(@cand_2,0);
  $this->choose_rep_gamfl(@cand_3,1);
  $this->choose_rep_gamfl(@cand_4,1);

  $this->clean_munu_4();
  return;
}

sub set_texf_4 {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};

  for ( my $i = 0 ; $i <= $#G ; $i++ ) {
    $G[$i] = &tex_gamma($G[$i]);
  }

  $this->{texf} = "F_4(x,y;$G[0],$G[1],$G[2],$G[3])^\{$F[0]$F[1]$F[2]$F[3]\}";
  return;
}

sub set_name_4 {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  $this->{name} = $this->{type} . "$G[0]$F[0]$G[1]$F[1]$G[2]$F[2]$G[3]$F[3]";
  return;
}

##      ##
#  F_4p  #
##      ##

sub get_cand_star_4p {
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-y) Γ4 Sf4(y-x) ]
  #  = Tr[ Γ2★ Sf1(x-x) Γ1★ Sf4(x-y) Γ4★ Sf3(y-y) Γ3★ Sf2(y-x) ]
  # This sub translates:
  #    (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4) -> (Γ2★,Γ1★,Γ4★,Γ3★)(F1,F4,F3,F2)
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  @G = &get_g5Hconjs(@G);
  @G = ($G[1],$G[0],$G[3],$G[2]);
  @F = ($F[0],$F[3],$F[2],$F[1]);
  return([@G],[@F]);
}

sub clean_munu_4p {
  my $this = shift;
  my @G = @{$this->{gammas}};
  $G[0] =~ s/mu\d/mu/;
  $G[1] =~ s/mu\d/mu/;
  $G[2] =~ s/mu\d/nu/;
  $G[3] =~ s/mu\d/nu/;
  $this->{gammas} = [@G];
  return;
}

sub set_rep_gamfl_4p {
  my $this = shift;
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-y) Γ4 Sf4(y-x) ]
  #   Cand:  (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4)
  #   Cand*: (Γ2★,Γ1★,Γ4★,Γ3★)(F1,F4,F3,F2)
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  my @cand_1 = ([@G],[@F]);
  my @cand_2 = &get_cand_star_4p(@cand_1);

  # Update to representative gamfl with conjf
  $this->choose_rep_gamfl(@cand_2,1);

  $this->clean_munu_4p();
  return;
}

sub set_texf_4p {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};

  for ( my $i = 0 ; $i <= $#G ; $i++ ) {
    $G[$i] = &tex_gamma($G[$i]);
  }

  $this->{texf} = "F_4'(x,y;$G[0],$G[1],$G[2],$G[3])^\{$F[0]$F[1]$F[2]$F[3]\}";
  return;
}

sub set_name_4p {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  $this->{name} = $this->{type} . "$G[0]$F[0]$G[1]$F[1]$G[2]$F[2]$G[3]$F[3]";
  return;
}

##      ##
#  F_3x  #
##      ##

sub get_cand_star_3x {
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-x) ]
  #  = Tr[ Γ2★ Sf1(x-x) Γ1★ Sf3(x-y) Γ3★ Sf2(y-x) ]
  # This sub translates:
  #    (Γ1,Γ2,Γ3)(F1,F2,F3) -> (Γ2★,Γ1★,Γ3★)(F1,F3,F2)
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  @G = &get_g5Hconjs(@G);
  @G = ($G[1],$G[0],$G[2]);
  @F = ($F[0],$F[2],$F[1]);
  return([@G],[@F]);
}

sub clean_munu_3x {
  my $this = shift;
  my @G = @{$this->{gammas}};
  $G[0] =~ s/mu\d/mu/;
  $G[1] =~ s/mu\d/mu/;
  $G[2] =~ s/mu\d/nu/;
  $this->{gammas} = [@G];
  return;
}

sub set_rep_gamfl_3x {
  my $this = shift;
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-y) Γ4 Sf4(y-x) ]
  #   Cand:  (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4)
  #   Cand*: (Γ2★,Γ1★,Γ3★)(F1,F3,F2)
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  my @cand_1 = ([@G],[@F]);
  my @cand_2 = &get_cand_star_3x(@cand_1);

  # Update to representative gamfl with conjf
  $this->choose_rep_gamfl(@cand_2,1);

  $this->clean_munu_3x();
  return;
}

sub set_texf_3x {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};

  for ( my $i = 0 ; $i <= $#G ; $i++ ) {
    $G[$i] = &tex_gamma($G[$i]);
  }

  $this->{texf} = "F_3(x,y;$G[0],$G[1],$G[2])^\{$F[0]$F[1]$F[2]\}";
  return;
}

sub set_name_3x {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  $this->{name} = $this->{type} . "$G[0]$F[0]$G[1]$F[1]$G[2]$F[2]";
  return;
}

##      ##
#  F_3y  #
##      ##

sub get_cand_star_3y {
  # Tr[ Γ1 Sf1(y-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) ]
  #  = Tr[ Γ2★ Sf1(y-y) Γ1★ Sf3(y-x) Γ3★ Sf2(x-y) ]
  # This sub translates:
  #    (Γ1,Γ2,Γ3)(F1,F2,F3) -> (Γ2★,Γ1★,Γ3★)(F1,F3,F2)
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  @G = &get_g5Hconjs(@G);
  @G = ($G[1],$G[0],$G[2]);
  @F = ($F[0],$F[2],$F[1]);
  return([@G],[@F]);
}

sub clean_munu_3y {
  my $this = shift;
  my @G = @{$this->{gammas}};
  $G[0] =~ s/mu\d/nu/;
  $G[1] =~ s/mu\d/nu/;
  $G[2] =~ s/mu\d/mu/;
  $this->{gammas} = [@G];
  return;
}

sub set_rep_gamfl_3y {
  my $this = shift;
  # Tr[ Γ1 Sf1(y-y) Γ2 Sf2(y-x) Γ3 Sf3(x-x) Γ4 Sf4(x-y) ]
  #   Cand:  (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4)
  #   Cand*: (Γ2★,Γ1★,Γ3★)(F1,F3,F2)
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  my @cand_1 = ([@G],[@F]);
  my @cand_2 = &get_cand_star_3y(@cand_1);

  # Update to representative gamfl with conjf
  $this->choose_rep_gamfl(@cand_2,1);

  $this->clean_munu_3y();
  return;
}

sub set_texf_3y {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};

  for ( my $i = 0 ; $i <= $#G ; $i++ ) {
    $G[$i] = &tex_gamma($G[$i]);
  }

  $this->{texf} = "F_3(y,x;$G[0],$G[1],$G[2])^\{$F[0]$F[1]$F[2]\}";
  return;
}

sub set_name_3y {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  $this->{name} = $this->{type} . "$G[0]$F[0]$G[1]$F[1]$G[2]$F[2]";
  return;
}

##     ##
#  F_2  #
##     ##

sub get_cand_star_2 {
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ]
  #  = Tr[ Γ1★ Sf1(x-y) Γ2★ Sf2(y-x) ]
  # This sub translates:
  #    (Γ1,Γ2)(F1,F2) -> (Γ1★,Γ2★)(F2,F1)
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  @G = &get_g5Hconjs(@G);
  @G = ($G[0],$G[1]);
  @F = ($F[1],$F[0]);
  return([@G],[@F]);
}

sub clean_munu_2 {
  my $this = shift;
  my @G = @{$this->{gammas}};
  $G[0] =~ s/mu\d/mu/;
  $G[1] =~ s/mu\d/nu/;
  $this->{gammas} = [@G];
  return;
}

sub set_rep_gamfl_2 {
  my $this = shift;
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ]
  #   Cand:  (Γ1,Γ2)(F1,F2)
  #   Cand*: (Γ1★,Γ2★)(F2,F1)
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  my @cand_1 = ([@G],[@F]);
  my @cand_2 = &get_cand_star_2(@cand_1);

  # Update to representative gamfl with conjf
  $this->choose_rep_gamfl(@cand_2,1);

  $this->clean_munu_2();
  return;
}

sub set_texf_2 {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};

  for ( my $i = 0 ; $i <= $#G ; $i++ ) {
    $G[$i] = &tex_gamma($G[$i]);
  }

#  $this->{texf} = "F_2(y,x;$G[0],$G[1])^\{$F[0]$F[1]\}";
  $this->{texf} = "F_2(x,y;$G[0],$G[1])^\{$F[0]$F[1]\}";
  return;
}

sub set_name_2 {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  $this->{name} = $this->{type} . "$G[0]$F[0]$G[1]$F[1]";
  return;
}

##      ##
#  F_1x  #
##      ##

sub get_cand_star_1x {
  # Tr[ Γ1 Sf1(x-x) ]
  #  = Tr[ Γ1★ Sf1(x-x) ]
  # This sub translates:
  #    (Γ1)(F1) -> (Γ1★)(F1)
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  @G = &get_g5Hconjs(@G);
  return([@G],[@F]);
}

sub clean_munu_1x {
  my $this = shift;
  my @G = @{$this->{gammas}};
  $G[0] =~ s/mu\d/mu/;
  $this->{gammas} = [@G];
  return;
}

sub set_rep_gamfl_1x {
  my $this = shift;
  # Tr[ Γ1 Sf1(x-y) ]
  #   Cand:  (Γ1)(F1)
  #   Cand*: (Γ1★)(F1)
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  my @cand_1 = ([@G],[@F]);
  my @cand_2 = &get_cand_star_1x(@cand_1);

  # Update to representative gamfl with conjf
  $this->choose_rep_gamfl(@cand_2,1);

  $this->clean_munu_1x();
  return;
}

sub set_texf_1x {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};

  for ( my $i = 0 ; $i <= $#G ; $i++ ) {
    $G[$i] = &tex_gamma($G[$i]);
  }

  $this->{texf} = "F_1(x;$G[0])^$F[0]";
  return;
}

sub set_name_1x {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  $this->{name} = $this->{type} . "$G[0]$F[0]";
  return;
}

##      ##
#  F_1y  #
##      ##

sub get_cand_star_1y {
  # Tr[ Γ1 Sf1(y-y) ]
  #  = Tr[ Γ1★ Sf1(y-y) ]
  # This sub translates:
  #    (Γ1)(F1) -> (Γ1★)(F1)
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  @G = &get_g5Hconjs(@G);
  return([@G],[@F]);
}

sub clean_munu_1y {
  my $this = shift;
  my @G = @{$this->{gammas}};
  $G[0] =~ s/mu\d/nu/;
  $this->{gammas} = [@G];
  return;
}

sub set_rep_gamfl_1y {
  my $this = shift;
  # Tr[ Γ1 Sf1(y-y) ]
  #   Cand:  (Γ1)(F1)
  #   Cand*: (Γ1★)(F1)
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  my @cand_1 = ([@G],[@F]);
  my @cand_2 = &get_cand_star_1y(@cand_1);

  # Update to representative gamfl with conjf
  $this->choose_rep_gamfl(@cand_2,1);

  $this->clean_munu_1y();
  return;
}

sub set_texf_1y {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};

  for ( my $i = 0 ; $i <= $#G ; $i++ ) {
    $G[$i] = &tex_gamma($G[$i]);
  }

  $this->{texf} = "F_1(y;$G[0])^$F[0]";
  return;
}

sub set_name_1y {
  my $this = shift;
  my @G = @{$this->{gammas}};
  my @F = @{$this->{flavors}};
  $this->{name} = $this->{type} . "$G[0]$F[0]";
  return;
}

#####################################

sub set_rep_gamfl {
  my $this = shift;
  return $this->set_rep_gamfl_4() if ( $this->{type} eq $type_4 );
  return $this->set_rep_gamfl_4p() if ( $this->{type} eq $type_4p );
  return $this->set_rep_gamfl_3x() if ( $this->{type} eq $type_3x );
  return $this->set_rep_gamfl_3y() if ( $this->{type} eq $type_3y );
  return $this->set_rep_gamfl_2() if ( $this->{type} eq $type_2 );
  return $this->set_rep_gamfl_1x() if ( $this->{type} eq $type_1x );
  return $this->set_rep_gamfl_1y() if ( $this->{type} eq $type_1y );
  die;
}

sub set_texf {
  my $this = shift;
  $this->set_texf_4() if ( $this->{type} eq $type_4 );
  $this->set_texf_4p() if ( $this->{type} eq $type_4p );
  $this->set_texf_3x() if ( $this->{type} eq $type_3x );
  $this->set_texf_3y() if ( $this->{type} eq $type_3y );
  $this->set_texf_2() if ( $this->{type} eq $type_2 );
  $this->set_texf_1x() if ( $this->{type} eq $type_1x );
  $this->set_texf_1y() if ( $this->{type} eq $type_1y );
  die unless defined( $this->{texf} );;
  if ( $this->{conjf} == 1 ) {
    $this->{texf} = "\{" . $this->{texf} . "\}^*";
  }
  return;
}

sub set_name {
  my $this = shift;
  $this->set_name_4() if ( $this->{type} eq $type_4 );
  $this->set_name_4p() if ( $this->{type} eq $type_4p );
  $this->set_name_3x() if ( $this->{type} eq $type_3x );
  $this->set_name_3y() if ( $this->{type} eq $type_3y );
  $this->set_name_2() if ( $this->{type} eq $type_2 );
  $this->set_name_1x() if ( $this->{type} eq $type_1x );
  $this->set_name_1y() if ( $this->{type} eq $type_1y );
  die unless defined( $this->{name} );;
  if ( $this->{conjf} == 1 ) {
#    $this->{name} .= "C";
  }
  return;
}

sub get_itype {
  my $this = shift;
  return 0 if ( $this->{type} eq $type_4 );
  return 1 if ( $this->{type} eq $type_4p );
  return 2 if ( $this->{type} eq $type_3x );
  return 3 if ( $this->{type} eq $type_3y );
  return 4 if ( $this->{type} eq $type_2 );
  return 5 if ( $this->{type} eq $type_1x );
  return 6 if ( $this->{type} eq $type_1y );
}

sub remove_conj_sign_flags {
  my $this = shift;
  $this->{conjf} = 0;
  $this->{sign} = 1;
  return;
}

sub get_n_oddmu {
  my $this = shift;
  my $mu = shift;
  my $nmu = 0;
  foreach my $G (@{$this->{gammas}}) {
    $nmu++ if ( $G =~ /$mu/ );
  }
  return ( $nmu % 2 );
}

sub copy_general {
  my $class = shift;
  my $this = shift;
  my $new = {
               type => $this->{type},
               gammas => $this->{gammas},
               flavors => $this->{flavors},
               texf => $this->{texf},
               name => $this->{name},
               conjf => $this->{conjf}
#               sign => undef
  };
  bless $new, $class;
  return $new;
}

1;
