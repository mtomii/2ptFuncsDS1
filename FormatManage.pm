package FormatManage;
use strict;
use warnings;
use RepManage;
use TraceManage;

our $light_strange = 1;# 0 -> physical, 1 -> SU(3) limit?
our $write_flag = &TraceManage::get_write_flag();

our $rtype_4    = &RepManage::get_type_4();
our $rtype_22   = &RepManage::get_type_22();
our $rtype_4p   = &RepManage::get_type_4p();
our $rtype_31xy = &RepManage::get_type_31xy();
our $rtype_31yx = &RepManage::get_type_31yx();
our $rtype_211  = &RepManage::get_type_211();
our $rtype_3x   = &RepManage::get_type_3x();
our $rtype_3y   = &RepManage::get_type_3y();
our $rtype_21x  = &RepManage::get_type_21x();
our $rtype_21y  = &RepManage::get_type_21y();
our $rtype_2    = &RepManage::get_type_2();

our $type_4   = "\$F_4\$";
our $type_22  = "\$F_{2,2}\$";
our $type_4p  = "\$F_4'\$";
our $type_31  = "\$F_{3,1}\$";
our $type_211 = "\$F_{2,1,1}\$";
our $type_3   = "\$F_3\$";
our $type_21  = "\$F_{2,1}\$";
our $type_2   = "\$F_2\$";

sub get_type_4 { return $type_4 }
sub get_type_22 { return $type_22 }
sub get_type_4p { return $type_4p }
sub get_type_31 { return $type_31 }
sub get_type_211 { return $type_211 }
sub get_type_3 { return $type_3 }
sub get_type_21 { return $type_21 }
sub get_type_2 { return $type_2 }

sub new {
  my $class = shift;
  my $tracef = shift;
  print "# FormatManage->new with $tracef\n" if ( $write_flag > 1 );
  my $this = {
               objs => [],
               nobjs => undef,
               nl => 0,
               ns => 0,
               nc => 0,
               type => undef,
               coef => undef
  };
  bless $this, $class;

  $this->set_type($tracef);# type set here
  my @gamfl = $this->get_gammas_flavors($tracef);
  foreach my $f (@{$gamfl[1]}) {
    my $ref = "n" . $f;
    $this->{$ref}++;
  }
  if ( $write_flag == 4 ) {
    print "# nl = $this->{nl}, ns = $this->{ns}, nc = $this->{nc}\n";
  }

  # set Rep1(, Rep2), object(s) of RepManage, in $this->{obj} = [Rep1(,Rep2)]
  # Each Rep is a representative of the same contractions
  # Rep1 and Rep2 are the same in the infinite statistic limit
  $this->set_rep_objs(@gamfl);
  my @objs = @{$this->{objs}};
  $this->{nobjs} = $#objs+1;
  if ( $write_flag > 1 ) {
    print "# Representatives:\n";
    foreach my $obj (@objs) {
      print "   $obj->{name}\n";
    }
  }

  return $this;
}

sub set_type {
  my $this = shift;
  my $tracef = shift;

  my $flag = 0;
  if ( $tracef =~ /Tr.+Tr.+Tr/ ) {# F_{2,1,1}
    $this->{type} = $type_211;
    print "# format type: $type_211\n" if ( $write_flag > 1 );
  } elsif ( $tracef =~ /Tr\[\s+(.+)\s+\]\s*Tr\[\s+(.+)\s+\]/ ) {
    # case of F_{2,2} || F_{3,1} || F_{2,1}
    my @TR1 = split('\s+',$1);
    my @TR2 = split('\s+',$2);
    if ( $#TR1 < $#TR2 ) {
      my @tmp = @TR1;
      @TR1 = @TR2;
      @TR2 = @tmp;
    }
    if ( $#TR1 * $#TR2 == 9 ) {# F_{2,2}
      $this->{type} = $type_22;
      print "# format type: $type_22\n" if ( $write_flag > 1 );
    } elsif ( $#TR1 * $#TR2 == 3 ) {# F_{2,1}
      $this->{type} = $type_21;
      print "# format type: $type_21\n" if ( $write_flag > 1 );
    } else {# F_{3,1}
      unless ( $#TR1 == 5 ) {
        print "@TR1\n";
        die;
      }
      $this->{type} = $type_31;
      print "# format type: $type_31\n" if ( $write_flag > 1 );
    }
  } elsif ( $tracef =~ /Tr\[\s+(.+)\s+\]/ ) {
    # case of F_4 || F_4' || F_3 || F_2
    my $tr = $1;
    my @TR = split('\s+',$tr);
    if ( $#TR == 3 ) {# F_2
      $this->{type} = $type_2;
      print "# format type: $type_2\n" if ( $write_flag > 1 );
    } elsif ( $#TR == 5 ) {# F_3
      $this->{type} = $type_3;
      print "# format type: $type_3\n" if ( $write_flag > 1 );
    } else {# F_4 || F_4'
      die unless ( $#TR == 7 );
      if ( $tr =~ /x-x/ ) {
        $this->{type} = $type_4p;
        print "# format type: $type_4p\n" if ( $write_flag > 1 );
      } else {
        $this->{type} = $type_4;
        print "# format type: $type_4\n" if ( $write_flag > 1 );
      }
    }
  }
  print "#---------------------------\n" if ( $write_flag > 1 );
  return;
}

sub flavor_clearner {
  my @flavors = @_;
  for ( my $i = 0 ; $i <= $#flavors ; $i++ ) {
    $flavors[$i] = $1 if ( $flavors[$i] =~ /S([udsc])/ );
    $flavors[$i] =~ s/[ud]/l/;
    $flavors[$i] =~ s/s/l/ if ( $light_strange > 0 );
  }
  return @flavors;
}

sub order_objs {
  my $rep1 = shift;
  my $rep2 = shift;
  my $name1 = $rep1->{name};
  my $name2 = $rep2->{name};
  return ($rep1) if ( $name1 eq $name2 );
  if ( $name1 =~ /S$/ ) {
    $name1 =~ s/S$//;
    return ($rep2) if ( $name2 eq $name1 );
  } elsif ( $name2 =~ /S$/ ) {
    $name2 =~ s/S$//;
    return ($rep1) if ( $name1 eq $name2 );
  }
  my $test1 = $name1;
  my $test2 = $name2;
  $test1 =~ s/C//g;
  $test2 =~ s/C//g;
  if ( $test1 eq $test2 ) {
    my @names1 = split('_',$name1);
    my @names2 = split('_',$name2);
    for ( my $i = 0 ; $i <= $#names1 ; $i++ ) {
      next if ( $names1[$i] eq $names2[$i] );
#      print "GGG $rep1->{name} \\= $rep2->{name}\n";
      return ($rep1,$rep2) if ( $names2[$i] =~ /C$/ );
      return ($rep2,$rep1) if ( $names1[$i] =~ /C$/ );
    }
  }

  my @t1 = @{$rep1->{traces}};
  my @t2 = @{$rep2->{traces}};

  for ( my $ir = 0 ; $ir <= $#t1 ; $ir++ ) {
    my $tr1 = $t1[$ir];
    my $tr2 = $t2[$ir];

    my @G1 = @{$tr1->{gammas}};
    my @F1 = @{$tr1->{flavors}};
    my @G2 = @{$tr2->{gammas}};
    my @F2 = @{$tr2->{flavors}};
    for ( my $i = 0 ; $i <= $#F1 ; $i++ ) {
      next if ( $F1[$i] eq $F2[$i] );
      my $iq1 = &TraceManage::get_iq($F1[$i]);
      my $iq2 = &TraceManage::get_iq($F2[$i]);
      if ( $iq1 < $iq2 ) {
        return ($rep1,$rep2);
      } else {
        return ($rep2,$rep1);
      }
    }
    for ( my $i = 0 ; $i <= $#G1 ; $i++ ) {
      next if ( $G1[$i] eq $G2[$i] );
      my $ig1 = &TraceManage::get_ig($G1[$i]);
      my $ig2 = &TraceManage::get_ig($G2[$i]);
      if ( $ig1 < $ig2 ) {
        return ($rep1,$rep2);
      } else {
        return ($rep2,$rep1);
      }
    }
  }
  die;
}

sub get_gammas_flavors_4 {
  my $tracef = shift;
  my @gammas = ();
  my @flavors = ();
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) Γ4 Sf4(y-x) ]
  my @TR = ();
  @TR = split('\s+',$1) if ( $tracef =~ /Tr\[\s+(.+)\s+\]/ );
  if ( $TR[1] =~ /x-y/ ) {
    @gammas  = ($TR[0],$TR[2],$TR[4],$TR[6]);
    @flavors = ($TR[1],$TR[3],$TR[5],$TR[7]);
  } elsif ( $TR[3] =~ /x-y/ ) {
    @gammas  = ($TR[2],$TR[4],$TR[6],$TR[0]);
    @flavors = ($TR[3],$TR[5],$TR[7],$TR[1]);
  } else {
    print "$tracef\n";
    die;
  }
  @flavors = &flavor_clearner(@flavors);
  if ( $write_flag > 1 ) {
    print "# Initial gammas:  @gammas\n";
    print "# Initial flavors: @flavors\n";
  }
  return ([@gammas],[@flavors]);
}

sub set_rep_objs_4 {
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) Γ4 Sf4(y-x) ]
  #   Cand1:  (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4), (Γ3,Γ4,Γ1,Γ2)(F3,F4,F1,F2)
  #   Cand2:  (Γ2,Γ3,Γ4,Γ1)(F2,F3,F4,F1), ...
  print "# sub set_rep_objs_4\n" if ( $write_flag > 1 );
  my $rep1 = RepManage->new($rtype_4,[@G],[@F]);
#  my @objs = ($rep1);
  @G = ($G[1],$G[2],$G[3],$G[0]);
  @F = ($F[1],$F[2],$F[3],$F[0]);
  my $rep2 = RepManage->new($rtype_4,[@G],[@F]);
  my @objs = &order_objs($rep1,$rep2);
#  push (@objs,$rep2) unless ( $rep1->{name} eq $rep2->{name} );
  $this->{objs} = [@objs];
  return;
}

sub get_gammas_flavors_22 {
  my $in = shift;
  my $tracef = $in;
  my @gammas = ();
  my @flavors = ();
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(x-y) Γ4 Sf4(y-x) ]
  while( $tracef =~ /^[^T]+Tr\[\s+([^\]]+)\s+\]/ ) {
    my @TR = split('\s+',$1);
    if ( $TR[1] =~ /x-y/ ) {
      push(@gammas, ($TR[0],$TR[2]));
      push(@flavors,($TR[1],$TR[3]));
    } elsif ( $TR[3] =~ /x-y/ ) {
      push(@gammas, ($TR[2],$TR[0]));
      push(@flavors,($TR[3],$TR[1]));
    } else {
      print "$tracef\n";
      die;
    }
    $tracef =~ s/^[^\]]+\]//;
  }
  @flavors = &flavor_clearner(@flavors);
  if ( $write_flag > 1 ) {
    print "# Initial gammas:  @gammas\n";
    print "# Initial flavors: @flavors\n";
  }
  return ([@gammas],[@flavors]);
}

sub set_rep_objs_22 {
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(x-y) Γ4 Sf4(y-x) ]
  #   Cand1:  (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4)
  #   Cand2:  (Γ2,Γ1,Γ4,Γ3)(F2,F1,F4,F3)
  print "# sub set_rep_objs_22\n" if ( $write_flag > 1 );
  my $rep1 = RepManage->new($rtype_22,[@G],[@F]);
#  my @objs = ($rep1);
  @G = ($G[1],$G[0],$G[3],$G[2]);
  @F = ($F[1],$F[0],$F[3],$F[2]);
  my $rep2 = RepManage->new($rtype_22,[@G],[@F]);
  my @objs = &order_objs($rep1,$rep2);
#  push (@objs,$rep2) unless ( $rep1->{name} eq $rep2->{name} );
  $this->{objs} = [@objs];
  return;
}

sub get_gammas_flavors_4p {
  my $tracef = shift;
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-y) Γ4 Sf4(y-x) ]
  my @TR = ();
  @TR = split('\s+',$1) if ( $tracef =~ /Tr\[\s+(.+)\s+\]/ );
  until ( $TR[1] =~ /x-x/ ) {
    unshift(@TR,pop(@TR));
    unshift(@TR,pop(@TR));
  }
  my @gammas  = ($TR[0],$TR[2],$TR[4],$TR[6]);
  my @flavors = ($TR[1],$TR[3],$TR[5],$TR[7]);
  @flavors = &flavor_clearner(@flavors);
  if ( $write_flag > 1 ) {
    print "# Initial gammas:  @gammas\n";
    print "# Initial flavors: @flavors\n";
  }
  return ([@gammas],[@flavors]);
}

sub set_rep_objs_4p {
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-y) Γ4 Sf4(y-x) ]
  #   Cand1: (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4)
  #   Cand2: (Γ3,Γ4,Γ1,Γ2)(F3,F4,F1,F2)
  print "# sub set_rep_objs_4p\n" if ( $write_flag > 1 );
  my $rep1 = RepManage->new($rtype_4p,[@G],[@F]);
#  my @objs = ($rep1);
  @G = ($G[2],$G[3],$G[0],$G[1]);
  @F = ($F[2],$F[3],$F[0],$F[1]);
  my $rep2 = RepManage->new($rtype_4p,[@G],[@F]);
  my @objs = &order_objs($rep1,$rep2);
#  push (@objs,$rep2) unless ( $rep1->{name} eq $rep2->{name} );
  $this->{objs} = [@objs];
  return;
}

sub get_gammas_flavors_31 {
  my $tracef = shift;
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-x) ] Tr[ Γ4 Sf4(y-y) ] ||
  # Tr[ Γ1 Sf1(y-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) ] Tr[ Γ4 Sf4(x-x) ]
  my @TR1 = ();
  my @TR2 = ();
  if ( $tracef =~ /Tr\[\s+(.+)\s+\]\s*Tr\[\s+(.+)\s+\]/ ) {
    @TR1 = split('\s+',$1);
    @TR2 = split('\s+',$2);
  } else {
    print "$tracef\n";
    die;
  }
  if ( $#TR1 < $#TR2 ) {
    my @tmp = @TR1;
    @TR1 = @TR2;
    @TR2 = @tmp;
  }
  until ( $TR1[1] =~ /x-x/ || $TR1[1] =~ /y-y/ ) {
    unshift(@TR1,pop(@TR1));
    unshift(@TR1,pop(@TR1));
  }
  my @gammas  = ($TR1[0],$TR1[2],$TR1[4],$TR2[0]);
  my @flavors = ($TR1[1],$TR1[3],$TR1[5],$TR2[1]);
  @flavors = &flavor_clearner(@flavors);
  if ( $write_flag > 1 ) {
    print "# Initial gammas:  @gammas\n";
    print "# Initial flavors: @flavors\n";
  }
  return ([@gammas],[@flavors]);
}

sub set_rep_objs_31 {
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-x) ] Tr[ Γ4 Sf4(y-y) ]
  #   Cand1: x(Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4)
  #   Cand2: y(Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4)
  print "# sub set_rep_objs_31\n" if ( $write_flag > 1 );
  my $rep1 = RepManage->new($rtype_31xy,[@G],[@F]);
  my $rep2 = RepManage->new($rtype_31yx,[@G],[@F]);
  my @objs = ($rep1,$rep2);
  $this->{objs} = [@objs];
  return;
}

sub get_gammas_flavors_211 {
  my $tracef = shift;
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(x-x) ] Tr[ Γ4 Sf4(y-y) ]
  my @TR1 = ();
  my @TR2 = ();
  my @TR3 = ();
  if ( $tracef =~ /Tr\[\s+(.+)\s+\]\s*Tr\[\s+(.+)\s+\]\s*Tr\[\s+(.+)\s+\]/ ) {
    my @comps = ($1,$2,$3);
    foreach my $trc (@comps) {
      die if ( ( $trc =~ /x-y/ && $trc =~ /x-x/ ) ||
               ( $trc =~ /x-x/ && $trc =~ /y-y/ ) ||
               ( $trc =~ /y-y/ && $trc =~ /x-y/ ) );
      @TR1 = split('\s+',$trc) if ( $trc =~ /x-y/ );
      @TR2 = split('\s+',$trc) if ( $trc =~ /x-x/ );
      @TR3 = split('\s+',$trc) if ( $trc =~ /y-y/ );
    }
    die if ( ( ($#TR1+1) * ($#TR2+1) * ($#TR3+1) == 0 ) );
  } else {
    print "$tracef\n";
    die;
  }
  until ( $TR1[1] =~ /x-y/ ) {
    unshift(@TR1,pop(@TR1));
    unshift(@TR1,pop(@TR1));
  }
  my @gammas  = ($TR1[0],$TR1[2],$TR2[0],$TR3[0]);
  my @flavors = ($TR1[1],$TR1[3],$TR2[1],$TR3[1]);
  @flavors = &flavor_clearner(@flavors);
  if ( $write_flag > 1 ) {
    print "# Initial gammas:  @gammas\n";
    print "# Initial flavors: @flavors\n";
  }
  return ([@gammas],[@flavors]);
}

sub set_rep_objs_211 {
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(x-x) ] Tr[ Γ4 Sf4(y-y) ]
  #   Cand1:  (Γ1,Γ2,Γ3,Γ4)(F1,F2,F3,F4)
  #   Cand2:  (Γ2,Γ1,Γ4,Γ3)(F2,F1,F4,F3)
  print "# sub set_rep_objs_211\n" if ( $write_flag > 1 );
  my $rep1 = RepManage->new($rtype_211,[@G],[@F]);
#  my @objs = ($rep1);
  @G = ($G[1],$G[0],$G[3],$G[2]);
  @F = ($F[1],$F[0],$F[3],$F[2]);
  my $rep2 = RepManage->new($rtype_211,[@G],[@F]);
  my @objs = &order_objs($rep1,$rep2);
#  push (@objs,$rep2) unless ( $rep1->{name} eq $rep2->{name} );
  $this->{objs} = [@objs];
  return;
}

sub get_gammas_flavors_3 {
  my $tracef = shift;
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-x) ] ||
  # Tr[ Γ1 Sf1(y-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) ]
  my @TR = ();
  if ( $tracef =~ /Tr\[\s+(.+)\s+\]/ ) {
    @TR = split('\s+',$1);
  } else {
    print "$tracef\n";
    die;
  }
  until ( $TR[1] =~ /x-x/ || $TR[1] =~ /y-y/ ) {
    unshift(@TR,pop(@TR));
    unshift(@TR,pop(@TR));
  }
  my @gammas  = ($TR[0],$TR[2],$TR[4]);
  my @flavors = ($TR[1],$TR[3],$TR[5]);
  @flavors = &flavor_clearner(@flavors);
  if ( $write_flag > 1 ) {
    print "# Initial gammas:  @gammas\n";
    print "# Initial flavors: @flavors\n";
  }
  return ([@gammas],[@flavors]);
}

sub set_rep_objs_3 {
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-x) ]
  #   Cand1: x(Γ1,Γ2,Γ3)(F1,F2,F3)
  #   Cand2: y(Γ1,Γ2,Γ3)(F1,F2,F3)
  print "# sub set_rep_objs_3\n" if ( $write_flag > 1 );
  my $rep1 = RepManage->new($rtype_3x,[@G],[@F]);
  my $rep2 = RepManage->new($rtype_3y,[@G],[@F]);
  my @objs = ($rep1,$rep2);
  $this->{objs} = [@objs];
  return;
}

sub get_gammas_flavors_21 {
  my $tracef = shift;
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(y-y) ] ||
  # Tr[ Γ1 Sf1(y-x) Γ2 Sf2(x-y) ] Tr[ Γ3 Sf3(x-x) ]
  my @TR1 = ();
  my @TR2 = ();
  if ( $tracef =~ /Tr\[\s+(.+)\s+\]\s*Tr\[\s+(.+)\s+\]/ ) {
    @TR1 = split('\s+',$1);
    @TR2 = split('\s+',$2);
  } else {
    print "$tracef\n";
    die;
  }
  if ( $#TR1 < $#TR2 ) {
    my @tmp = @TR1;
    @TR1 = @TR2;
    @TR2 = @tmp;
  }
  my $x1 = "x-y";
  $x1 = "y-x" if ( $TR2[1] =~ /x-x/ );
  until ( $TR1[1] =~ /$x1/ ) {
    unshift(@TR1,pop(@TR1));
    unshift(@TR1,pop(@TR1));
  }
  my @gammas  = ($TR1[0],$TR1[2],$TR2[0]);
  my @flavors = ($TR1[1],$TR1[3],$TR2[1]);
  @flavors = &flavor_clearner(@flavors);
  if ( $write_flag > 1 ) {
    print "# Initial gammas:  @gammas\n";
    print "# Initial flavors: @flavors\n";
  }
  return ([@gammas],[@flavors]);
}

sub set_rep_objs_21 {
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(y-y) ]
  #   Cand1:  y(Γ1,Γ2,Γ3)(F1,F2,F3)
  #   Cand2:  x(Γ2,Γ1,Γ3)(F2,F1,F3)
  print "# sub set_rep_objs_21\n" if ( $write_flag > 1 );
  my $rep1 = RepManage->new($rtype_21y,[@G],[@F]);
  @G = ($G[1],$G[0],$G[2]);
  @F = ($F[1],$F[0],$F[2]);
  my $rep2 = RepManage->new($rtype_21x,[@G],[@F]);
  my @objs = ($rep2,$rep1);
  $this->{objs} = [@objs];
  return;
}

sub get_gammas_flavors_2 {
  my $tracef = shift;
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ]
  my @TR = ();
  if ( $tracef =~ /Tr\[\s+(.+)\s+\]/ ) {
    @TR = split('\s+',$1);
  } else {
    print "$tracef\n";
    die;
  }
  until ( $TR[1] =~ /x-y/ ) {
    unshift(@TR,pop(@TR));
    unshift(@TR,pop(@TR));
  }
  my @gammas  = ($TR[0],$TR[2]);
  my @flavors = ($TR[1],$TR[3]);
  @flavors = &flavor_clearner(@flavors);
  if ( $write_flag > 1 ) {
    print "# Initial gammas:  @gammas\n";
    print "# Initial flavors: @flavors\n";
  }
  return ([@gammas],[@flavors]);
}

sub set_rep_objs_2 {
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ]
  #   Cand1:  y(Γ1,Γ2)(F1,F2)
  #   Cand2:  x(Γ2,Γ1)(F2,F1)
  print "# sub set_rep_objs_2\n" if ( $write_flag > 1 );
  my $rep1 = RepManage->new($rtype_2,[@G],[@F]);
#  my @objs = ($rep1);
  @G = ($G[1],$G[0]);
  @F = ($F[1],$F[0]);
  my $rep2 = RepManage->new($rtype_2,[@G],[@F]);
  my @objs = &order_objs($rep1,$rep2);
#  push (@objs,$rep2) unless ( $rep1->{name} eq $rep2->{name} );
  $this->{objs} = [@objs];
  return;
}

sub get_gammas_flavors {
  my $this = shift;
  my $tracef = shift;
  return &get_gammas_flavors_4($tracef) if ( $this->{type} eq $type_4 );
  return &get_gammas_flavors_22($tracef) if ( $this->{type} eq $type_22 );
  return &get_gammas_flavors_4p($tracef) if ( $this->{type} eq $type_4p );
  return &get_gammas_flavors_31($tracef) if ( $this->{type} eq $type_31 );
  return &get_gammas_flavors_211($tracef) if ( $this->{type} eq $type_211 );
  return &get_gammas_flavors_3($tracef) if ( $this->{type} eq $type_3 );
  return &get_gammas_flavors_21($tracef) if ( $this->{type} eq $type_21 );
  return &get_gammas_flavors_2($tracef) if ( $this->{type} eq $type_2 );
}

sub set_rep_objs {
  my $this = shift;
  my @gamfl = @_;
  return $this->set_rep_objs_4(@gamfl) if ( $this->{type} eq $type_4 );
  return $this->set_rep_objs_22(@gamfl) if ( $this->{type} eq $type_22 );
  return $this->set_rep_objs_4p(@gamfl) if ( $this->{type} eq $type_4p );
  return $this->set_rep_objs_31(@gamfl) if ( $this->{type} eq $type_31 );
  return $this->set_rep_objs_211(@gamfl) if ( $this->{type} eq $type_211 );
  return $this->set_rep_objs_3(@gamfl) if ( $this->{type} eq $type_3 );
  return $this->set_rep_objs_21(@gamfl) if ( $this->{type} eq $type_21 );
  return $this->set_rep_objs_2(@gamfl) if ( $this->{type} eq $type_2 );
}

sub get_texf_rep {
  my $this = shift;
  my $rep = $this->{objs}->[0];
  return $rep->{texf};
}

sub get_name_rep {
  my $this = shift;
  my $rep = $this->{objs}->[0];
  print "get_name_rep: $rep->{name}\n";
  return $rep->{name};
}

sub get_sign_rep {
  my $this = shift;
  my $rep = $this->{objs}->[0];
  print "get_sign_rep: $rep->{sign}\n";
  return $rep->{sign};
}

sub get_itype {
  my $this = shift;
  return 0 if ( $this->{type} eq $type_4 );
  return 1 if ( $this->{type} eq $type_22 );
  return 2 if ( $this->{type} eq $type_4p );
  return 3 if ( $this->{type} eq $type_31 );
  return 4 if ( $this->{type} eq $type_211 );
  return 5 if ( $this->{type} eq $type_3 );
  return 6 if ( $this->{type} eq $type_21 );
  return 7 if ( $this->{type} eq $type_2 );
}

sub remove_conj_sign_flags {
  my $this = shift;
  my @reps = @{$this->{objs}};
  for ( my $i = 0 ; $i <= $#reps ; $i++ ) {
    my $rep = $reps[$i];
    $rep->remove_conj_sign_flags();
    $reps[$i] = $rep;
  }
  $this->{objs} = [@reps];
  return;
}

sub copy_general {
  my $class = shift;
  my $this = shift;
  my $new = {
               objs => [],
               nobjs => $this->{nobjs},
               nl => $this->{nl},
               ns => $this->{ns},
               nc => $this->{nc},
               type => $this->{type}
  };

  bless $new, $class;

  my @new_reps = ();
  foreach my $rep(@{$this->{objs}}) {
    my $new_rep = RepManage->copy_general($rep);
    push(@new_reps,$new_rep);
  }
  $new->{objs} = [@new_reps];
  return $new;
}

1;
