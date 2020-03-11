package RepManage;
use strict;
use warnings;
use TraceManage;

our $write_flag = &TraceManage::get_write_flag();

our $trtype_4   = &TraceManage::get_type_4();
our $trtype_4p  = &TraceManage::get_type_4p();
our $trtype_3x  = &TraceManage::get_type_3x();
our $trtype_3y  = &TraceManage::get_type_3y();
our $trtype_2   = &TraceManage::get_type_2();
our $trtype_1x  = &TraceManage::get_type_1x();
our $trtype_1y  = &TraceManage::get_type_1y();

our $Rep = "Rep";
our $type_4    = $Rep . "_" . $trtype_4;
our $type_22   = $Rep . "_" . $trtype_2  . "_" . $trtype_2;
our $type_4p   = $Rep . "_" . $trtype_4p;
our $type_31xy = $Rep . "_" . $trtype_3x . "_" . $trtype_1y;
our $type_31yx = $Rep . "_" . $trtype_3y . "_" . $trtype_1x;
our $type_211  = $Rep . "_" . $trtype_2  . "_" . $trtype_1x . "_" . $trtype_1y;
our $type_3x   = $Rep . "_" . $trtype_3x;
our $type_3y   = $Rep . "_" . $trtype_3y;
our $type_21x  = $Rep . "_" . $trtype_2  . "_" . $trtype_1x;
our $type_21y  = $Rep . "_" . $trtype_2  . "_" . $trtype_1y;
our $type_2    = $Rep . "_" . $trtype_2;
sub get_type_4 { return $type_4 }
sub get_type_22 { return $type_22 }
sub get_type_4p { return $type_4p }
sub get_type_31xy { return $type_31xy }
sub get_type_31yx { return $type_31yx }
sub get_type_211 { return $type_211 }
sub get_type_3x { return $type_3x }
sub get_type_3y { return $type_3y }
sub get_type_21x { return $type_21x }
sub get_type_21y { return $type_21y }
sub get_type_2 { return $type_2 }

sub new {
  my $class = shift;
  my $type = shift;
  my @gamfl = @_;
  if ( $write_flag > 1 ) {
    print "## RepManage->new with rep type $type\n";
    my @G = @{$gamfl[0]};
    print "    gammas:  @G\n";
    my @F = @{$gamfl[1]};
    print "    flavors: @F\n";
  }
  my $this = {
               type => $type,
               traces => undef,
               texf => undef,
               name => undef,
               sign => undef,
               conjf => 0
  };

  bless $this, $class;

  # List up exactly the same traces and choose one as a representative
  $this->set_rep(@gamfl);
  $this->{texf} = "\\Big\\langle " . $this->{texf} . "\\Big\\rangle";
  if ( $this->{conjf} == 1 ) {
    $this->{texf} .= "^*";
    #$this->{name} .= "S"; needed if imaginary part is also concerned but need more treatment in analysis code to take care of this symbol
  }
  $this->set_sign();
  if ( $write_flag > 1 ) {
    print "## texf of Rep: $this->{texf}\n";
    print "## name of Rep: $this->{name}\n";
    print "## sign: $this->{sign}\n";
  }

  return $this;
}

sub set_rep_4 {
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) Γ4 Sf4(y-x) ]
  my $this = shift;
  my @gamfl = @_;

  my $tr4 = TraceManage->new($trtype_4,@gamfl);

  $this->{traces} = [$tr4];
  $this->nonCC_first_tr();
  $this->{texf} = $tr4->{texf};
  $this->{name} = $tr4->{name};
  $this->{name} .= "C" if ( $tr4->{conjf} == 1 );
  return $this;
}

sub order_rep_22 {
  my $tr2A = shift;
  my $tr2B = shift;
  my @GA = @{$tr2A->{gammas}};
  my @FA = @{$tr2A->{flavors}};
  my @GB = @{$tr2B->{gammas}};
  my @FB = @{$tr2B->{flavors}};
  for ( my $i = 0 ; $i <= 1 ; $i++ ) {
    next if ( $FA[$i] eq $FB[$i] );
    my $iqA = &TraceManage::get_iq($FA[$i]);
    my $iqB = &TraceManage::get_iq($FB[$i]);
    return ($tr2A,$tr2B) if ( $iqA < $iqB );
    return ($tr2B,$tr2A)
  }
  for ( my $i = 0 ; $i <= 1 ; $i++ ) {
    next if ( $GA[$i] eq $GB[$i] );
    my $igA = &TraceManage::get_ig($GA[$i]);
    my $igB = &TraceManage::get_ig($GB[$i]);
    return ($tr2A,$tr2B) if ( $igA < $igB );
    return ($tr2B,$tr2A)
  }
  return ($tr2A,$tr2B);
}

sub set_rep_22 {
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] [ Γ3 Sf3(x-y) Γ4 Sf4(y-x) ]
  my $this = shift;
  my @gamfl = @_;
  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  my @GA = ($G[0],$G[1]);
  my @FA = ($F[0],$F[1]);
  my @GB = ($G[2],$G[3]);
  my @FB = ($F[2],$F[3]);

  my $tr2A = TraceManage->new($trtype_2,[@GA],[@FA]);
  my $tr2B = TraceManage->new($trtype_2,[@GB],[@FB]);
  my @trs = &order_rep_22($tr2A,$tr2B);
  $this->{traces} = [@trs];
  $this->nonCC_first_tr();# $this->{traces} might be updated
  @trs = @{$this->{traces}};
  $tr2A = $trs[0];
  $tr2B = $trs[1];
  my $csignA = "";
  my $csignB = "";
  $csignA = "C" if ( $tr2A->{conjf} == 1 );
  $csignB = "C" if ( $tr2B->{conjf} == 1 );
  $this->{texf} = $tr2A->{texf} . " \\cdot " . $tr2B->{texf};
  $this->{name} = $tr2A->{name} . $csignA . "_" . $tr2B->{name} . $csignB;
  return $this;
}

sub set_rep_4p {
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-y) Γ4 Sf4(y-x) ]
  my $this = shift;
  my @gamfl = @_;

  my $tr4p = TraceManage->new($trtype_4p,@gamfl);

  $this->{traces} = [$tr4p];
  $this->nonCC_first_tr();
  $this->{texf} = $tr4p->{texf};
  $this->{name} = $tr4p->{name};
  $this->{name} .= "C" if ( $tr4p->{conjf} == 1 );
  return $this;
}

sub set_rep_31xy {
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-x) ] Tr[ Γ4 Sf4(y-y) ]
  my $this = shift;
  my @gamfl = @_;

  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  my @G3 = ($G[0],$G[1],$G[2]);
  my @F3 = ($F[0],$F[1],$F[2]);
  my @G1 = ($G[3]);
  my @F1 = ($F[3]);
  my $tr3 = TraceManage->new($trtype_3x,[@G3],[@F3]);
  my $tr1 = TraceManage->new($trtype_1y,[@G1],[@F1]);
  my $csign3 = "";
  my $csign1 = "";
  $csign3 = "C" if ( $tr3->{conjf} == 1 );
  $csign1 = "C" if ( $tr1->{conjf} == 1 );

  $this->{traces} = [$tr3,$tr1];
  $this->nonCC_first_tr();
  $this->{texf} = $tr3->{texf} . " \\cdot " . $tr1->{texf};
  $this->{name} = $tr3->{name} . $csign3 . "_" . $tr1->{name} . $csign1;
  return $this;
}

sub set_rep_31yx {
  # Tr[ Γ1 Sf1(y-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) ] Tr[ Γ4 Sf4(x-x) ]
  my $this = shift;
  my @gamfl = @_;

  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  my @G3 = ($G[0],$G[1],$G[2]);
  my @F3 = ($F[0],$F[1],$F[2]);
  my @G1 = ($G[3]);
  my @F1 = ($F[3]);
  my $tr3 = TraceManage->new($trtype_3y,[@G3],[@F3]);
  my $tr1 = TraceManage->new($trtype_1x,[@G1],[@F1]);
  my $csign3 = "";
  my $csign1 = "";
  $csign3 = "C" if ( $tr3->{conjf} == 1 );
  $csign1 = "C" if ( $tr1->{conjf} == 1 );

  $this->{traces} = [$tr3,$tr1];
  $this->nonCC_first_tr();
  $this->{texf} = $tr3->{texf} . " \\cdot " . $tr1->{texf};
  $this->{name} = $tr3->{name} . $csign3 . "_" . $tr1->{name} . $csign1;
  return $this;
}

sub set_rep_211 {
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(x-x) ] Tr[ Γ4 Sf4(y-y) ]
  my $this = shift;
  my @gamfl = @_;

  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  my @G2 = ($G[0],$G[1]);
  my @F2 = ($F[0],$F[1]);
  my @G1x = ($G[2]);
  my @F1x = ($F[2]);
  my @G1y = ($G[3]);
  my @F1y = ($F[3]);
  my $tr2  = TraceManage->new($trtype_2,[@G2],[@F2]);
  my $tr1x = TraceManage->new($trtype_1x,[@G1x],[@F1x]);
  my $tr1y = TraceManage->new($trtype_1y,[@G1y],[@F1y]);
  my $csign2 = "";
  my $csign1x = "";
  my $csign1y = "";
  $csign2 = "C" if ( $tr2->{conjf} == 1 );
  $csign1x = "C" if ( $tr1x->{conjf} == 1 );
  $csign1y = "C" if ( $tr1y->{conjf} == 1 );

  $this->{traces} = [$tr2,$tr1x,$tr1y];
  $this->nonCC_first_tr();
  $this->{texf} = $tr2->{texf} . " \\cdot " . $tr1x->{texf}
                               . " \\cdot " . $tr1y->{texf};
  $this->{name} = $tr2->{name} . $csign2 . "_" . $tr1x->{name} . $csign1x
                . "_" . $tr1y->{name} . $csign1y;
  return $this;
}

sub set_rep_3x {
  # Tr[ Γ1 Sf1(x-x) Γ2 Sf2(x-y) Γ3 Sf3(y-x) ]
  my $this = shift;
  my @gamfl = @_;

  my $tr3 = TraceManage->new($trtype_3x,@gamfl);

  $this->{traces} = [$tr3];
  $this->nonCC_first_tr();
  $this->{texf} = $tr3->{texf};
  $this->{name} = $tr3->{name};
  $this->{name} .= "C" if ( $tr3->{conjf} == 1 );
  return $this;
}

sub set_rep_3y {
  # Tr[ Γ1 Sf1(y-y) Γ2 Sf2(y-x) Γ3 Sf3(x-y) ]
  my $this = shift;
  my @gamfl = @_;

  my $tr3 = TraceManage->new($trtype_3y,@gamfl);

  $this->{traces} = [$tr3];
  $this->nonCC_first_tr();
  $this->{texf} = $tr3->{texf};
  $this->{name} = $tr3->{name};
  $this->{name} .= "C" if ( $tr3->{conjf} == 1 );
  return $this;
}

sub set_rep_21x {
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(x-x) ]
  my $this = shift;
  my @gamfl = @_;

  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  my @G2 = ($G[0],$G[1]);
  my @F2 = ($F[0],$F[1]);
  my @G1 = ($G[2]);
  my @F1 = ($F[2]);
  my $tr2 = TraceManage->new($trtype_2, [@G2],[@F2]);
  my $tr1 = TraceManage->new($trtype_1x,[@G1],[@F1]);
  my $csign2 = "";
  my $csign1 = "";
  $csign2 = "C" if ( $tr2->{conjf} == 1 );
  $csign1 = "C" if ( $tr1->{conjf} == 1 );

  $this->{traces} = [$tr2,$tr1];
  $this->nonCC_first_tr();
  $this->{texf} = $tr2->{texf} . " \\cdot " . $tr1->{texf};
  $this->{name} = $tr2->{name} . $csign2 . "_" . $tr1->{name} . $csign1;
  return $this;
}

sub set_rep_21y {
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ] Tr[ Γ3 Sf3(y-y) ]
  my $this = shift;
  my @gamfl = @_;

  my @G = @{$gamfl[0]};
  my @F = @{$gamfl[1]};
  my @G2 = ($G[0],$G[1]);
  my @F2 = ($F[0],$F[1]);
  my @G1 = ($G[2]);
  my @F1 = ($F[2]);
  my $tr2 = TraceManage->new($trtype_2, [@G2],[@F2]);
  my $tr1 = TraceManage->new($trtype_1y,[@G1],[@F1]);
  my $csign2 = "";
  my $csign1 = "";
  $csign2 = "C" if ( $tr2->{conjf} == 1 );
  $csign1 = "C" if ( $tr1->{conjf} == 1 );

  $this->{traces} = [$tr2,$tr1];
  $this->nonCC_first_tr();
  $this->{texf} = $tr2->{texf} . " \\cdot " . $tr1->{texf};
  $this->{name} = $tr2->{name} . $csign2 . "_" . $tr1->{name} . $csign1;
  return $this;
}

sub set_rep_2 {
  # Tr[ Γ1 Sf1(x-y) Γ2 Sf2(y-x) ]
  my $this = shift;
  my @gamfl = @_;

  my $tr2 = TraceManage->new($trtype_2,@gamfl);

  $this->{traces} = [$tr2];
  $this->nonCC_first_tr();
  $this->{texf} = $tr2->{texf};
  $this->{name} = $tr2->{name};
  $this->{name} .= "C" if ( $tr2->{conjf} == 1 );
  return $this;
}

sub set_rep {
  my $this = shift;
  my @gamfl = @_;
  return $this->set_rep_4(@gamfl) if ( $this->{type} eq $type_4 );
  return $this->set_rep_22(@gamfl) if ( $this->{type} eq $type_22 );
  return $this->set_rep_4p(@gamfl) if ( $this->{type} eq $type_4p );
  return $this->set_rep_31xy(@gamfl) if ( $this->{type} eq $type_31xy );
  return $this->set_rep_31yx(@gamfl) if ( $this->{type} eq $type_31yx );
  return $this->set_rep_211(@gamfl) if ( $this->{type} eq $type_211 );
  return $this->set_rep_3x(@gamfl) if ( $this->{type} eq $type_3x );
  return $this->set_rep_3y(@gamfl) if ( $this->{type} eq $type_3y );
  return $this->set_rep_21x(@gamfl) if ( $this->{type} eq $type_21x );
  return $this->set_rep_21y(@gamfl) if ( $this->{type} eq $type_21y );
  return $this->set_rep_2(@gamfl) if ( $this->{type} eq $type_2 );
  die;
}

sub set_sign {
  my $this = shift;
  $this->{sign} = 1;
  foreach my $tr (@{$this->{traces}}) {
    $this->{sign} *= $tr->{sign};
  }
  return;
}

sub nonCC_first_tr {
  my $this = shift;
  my $fsttr = $this->{traces}->[0];
  return if ( $fsttr->{conjf} == 0 );
  die unless ( $fsttr->{conjf} == 1 );
  my @traces = @{$this->{traces}};
  for ( my $i = 0 ; $i <= $#traces ; $i++ ) {
    my $tr = $traces[$i];
    if ( $tr->{conjf} == 0 ) {
      $tr->{conjf} = 1;
      $tr->{texf} = "\{" . $tr->{texf} . "\}^*";
    } elsif ( $tr->{conjf} == 1 ) {
      $tr->{conjf} = 0;
      die unless ( $tr->{texf} =~ /^\s*\{.+\}\^\*$/ );
      $tr->{texf} =~ s/^\{(.+)\}\^\*$/$1/;
    }
    $traces[$i] = $tr;
  }
  $this->{traces} = [@traces];
  $this->{conjf} = 1;
}

sub remove_conj_sign_flags {
  my $this = shift;
  $this->{name} =~ s/S$//;
  $this->{sign} = 1;
  $this->{conjf} = 0;
  $this->{texf} =~ s/\^\*$//;
  return;
}

sub copy_general {
  my $class = shift;
  my $this = shift;
  my $new = {
               type => $this->{type},
               traces => [],
               texf => undef,
               name => undef
#               sign => undef,
#               conjf => 0
  };
  bless $new, $class;

  my $texf = $this->{texf};
  $texf =~ s/\^\*$//;
  $new->{texf} = $texf;
  my $name = $this->{name};
  $name =~ s/S$//;
  $new->{name} = $name;

  my @new_trs = ();
  foreach my $tr(@{$this->{traces}}) {
    my $new_tr = TraceManage->copy_general($tr);
    push (@new_trs,$new_tr);
  }
  $new->{traces} = [@new_trs];
  return $new;
}

1;
