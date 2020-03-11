package CorrManage;
use strict;
use warnings;
use Functions;

our @quarks = ("u","d","s","c");
our @fields = ("ub","db","sb","cb","uv","dv","sv","cv");
our @positions = ("x","y","z","w");

sub new {
  my $class = shift;
  my @Oprs = @_;

  my $this = {
               corr => undef,
               coef => [1,0],
               gamma => [],
               num_f => undef,
               tracef => undef,
               texf => undef
              };
  bless $this, $class;

  die "More than 4pt funcs are not covered\n" if ($#Oprs>3);
  for ( my $io = 0 ; $io <= $#Oprs ; $io++ ) {
    $this->{corr} .= &get_func_O($Oprs[$io],$positions[$io]);
  }

  while ( $this->{corr} =~ /\s+(\d+)\s+/ ) {
    $this->{coef}->[0] *= $1;
    $this->{coef}->[1] *= $1;
    $this->{corr} =~ s/\s+$1\s+/ /;
  }
  while ( $this->{corr} =~ /\s-/ ) {
    $this->{corr} =~ s/\s-//;
    $this->{coef}->[0] *= -1;
    $this->{coef}->[1] *= -1;
  }
  while ( $this->{corr} =~ /\si/ ) {
    $this->{corr} =~ s/\si//;
    my $tmp = $this->{coef}->[0];
    $this->{coef}->[0] = - $this->{coef}->[1];
    $this->{coef}->[1] = $tmp;
  }
  while ( $this->{corr} =~ /\s(g[^\^]+\^\d+\^\d+)/ ) {
    my @tmp = @{$this->{gamma}};
    push(@tmp,$1);
    $this->{gamma} = [@tmp];
    $this->{corr} =~ s/\sg[^\^]+\^\d+\^\d+//;
  }

  my $test = $this->{corr};
  $test =~ s/^\s+//;
  my @testa = (split('\s+',$test));
  $this->{num_f} = $#testa + 1;
  return $this;
}


sub get_func_O {
  my $O = shift;
  my $x = shift;
  my $FO = $O;
  foreach my $f (@fields) {
    $FO =~ s/\s$f/ \/$f$x/g;
  }
  return $FO;
}

sub signManage {
  my @cont = @_;
  my @result = ();
  foreach my $con(@cont) {
    my $sign = 1;
    while ( $con =~ /-\s/ ) {
      $con =~ s/-\s//;
      $sign *= -1;
    }
    $con =~ s/\s\s+/ /g;
    push(@result," +$con") if ( $sign == 1 );
    push(@result," -$con") if ( $sign == -1 );
  }
  return @result;
}

sub trace_contraction {
  my $this = shift;
  my $type_label = shift;
  my $show_flag = ( shift || 0 );
  my @cont = ($this->{corr} . " E");
  &show_contract([@cont],"AAA") unless ( $show_flag == 0 );

  my $ctr = 0;

  my $flag = 1;
  while( $flag == 1 ) {
    my @tmp = ();
    foreach my $cnt (@cont) {
      if ( $cnt =~ /\/([@quarks])b([@positions])_(\d)\^(\d)([^\/]+)$/ ) {
        my $q = $1;
        my $x = $2;
        my $c = $3;
        my $s = $4;
        my $rem = $5;
        if ( $cnt =~ /\sE/ ) {
          my $next = $cnt;
          my $diff = "D$q" . "b$x" . "_$c^$s" . $rem;
          $next =~ s/\/[^\/]+$/$diff/;
          $next =~ s/\sE//;

          # (\del/\del \eta) \bar\eta D \eta brings negative sign
          $next = " -" . $next;

          # test if E is needed for next
          my $test = $next;
          my $ctrD = 0;
          my $ctre = 0;
          while ( $test =~ /\// ) {
            $test =~ s/\///;
            $ctre++;
          }
          while ( $test =~ /D/ ) {
            $test =~ s/D//;
            $ctrD++;
          }
          $next .= " E" if ( $ctre > $ctrD );
          
          push(@tmp,$next);
        }
        my $cand = "D$q" . "v";
        if ( $cnt =~ /$cand([@positions])_(\d)\^(\d).*$cand([@positions])_(\d)\^(\d)/ ) {
          my $y2 = $1;
          my $c2 = $2;
          my $s2 = $3;
          my $y3 = $4;
          my $c3 = $5;
          my $s3 = $6;

          # operate on former cand
          my $next = $cnt;
          my $test = $next;
          $test =~ s/$cand(.+$cand).+//;
          my $sign = 1;
          while ( $test =~ /D/ ) {
            $test =~ s/D//;
            $sign *= -1;
          }
          $next = " -" . $next if ( $sign == -1 );
          my $diff = "S$q($y2-$x)_$c2$c^$s2$s" . $rem;
          $next =~ s/\/[^\/]+$/$diff/;
          $next =~ s/($cand)[@positions]_\d\^\d(.*$cand)([@positions]_\d\^\d)/$2$3/;
          push(@tmp,$next);

          # operate on latter cand
          $next = $cnt;
          $test = $next;
          $test =~ s/($cand.+)$cand.+/$1/;
          $sign = 1;
          while ( $test =~ /D/ ) {
            $test =~ s/D//;
            $sign *= -1;
          }
          $next = " -" . $next if ( $sign == -1 );
          $diff = "S$q($y3-$x)_$c3$c^$s3$s" . $rem;
          $next =~ s/\/[^\/]+$/$diff/;
          $next =~ s/($cand)([@positions]_\d\^\d.*)($cand)[@positions]_\d\^\d/$1$2/;
          push(@tmp,$next);
        } elsif ( $cnt =~ /$cand([@positions])_(\d)\^(\d)/ ) {
          my $y = $1;
          my $c2 = $2;
          my $s2 = $3;
          my $next = $cnt;

          my $test = $next;
          my $sign = 1;
          $test =~ s/$cand(.+)$//;
          while ( $test =~ /D/ ) {
            $test =~ s/D//;
            $sign *= -1;
          }
          $next = " -" . $next if ( $sign == -1 );

          my $diff = "S$q($y-$x)_$c2$c^$s2$s" . $rem;
          $next =~ s/\/[^\/]+$/$diff/;
          $next =~ s/$cand([@positions])_\d\^\d//;

          push(@tmp,$next);
        }
      } elsif ( $cnt =~ /\/([@quarks])v([@positions])_(\d)\^(\d)([^\/]+)$/ ) {
        my $q = $1;
        my $x = $2;
        my $c = $3;
        my $s = $4;
        my $rem = $5;
        if ( $cnt =~ /\sE/ ) {
          my $next = $cnt;
          my $diff = "D$q" . "v$x" . "_$c^$s" . $rem;
#          $next =~ s/\/[@quarks]v[@positions]_\d\^\d[^\/]+$/$diff/;
          $next =~ s/\/[^\/]+$/$diff/;
          $next =~ s/\sE//;

          # test if E is needed for next
          my $test = $next;
          my $ctrD = 0;
          my $ctre = 0;
          while ( $test =~ /\// ) {
            $test =~ s/\///;
            $ctre++;
          }
          while ( $test =~ /D/ ) {
            $test =~ s/D//;
            $ctrD++;
          }
          $next .= " E" if ( $ctre > $ctrD );

          push(@tmp,$next);
        }
        my $cand = "D$q" . "b";
        if ( $cnt =~ /$cand([@positions])_(\d)\^(\d).*$cand([@positions])_(\d)\^(\d)/ ) {
          my $y2 = $1;
          my $c2 = $2;
          my $s2 = $3;
          my $y3 = $4;
          my $c3 = $5;
          my $s3 = $6;

          # operate on former cand
          my $next = $cnt;
          my $test = $next;
          $test =~ s/$cand(.+$cand).+//;
          my $sign = 1;
          while ( $test =~ /D/ ) {
            $test =~ s/D//;
            $sign *= -1;
          }
          $next = " -" . $next if ( $sign == -1 );
          my $diff = "S$q($x-$y2)_$c$c2^$s$s2" . $rem;
          $next =~ s/\/[^\/]+$/$diff/;
          $next =~ s/($cand)[@positions]_\d\^\d(.*$cand)([@positions]_\d\^\d)/$2$3/;
          push(@tmp,$next);

          # operate on latter cand
          $next = $cnt;
          $test = $next;
          $test =~ s/($cand.+)$cand.+/$1/;
          $sign = 1;
          while ( $test =~ /D/ ) {
            $test =~ s/D//;
            $sign *= -1;
          }
          $next = " -" . $next if ( $sign == -1 );
          $diff = "S$q($x-$y3)_$c$c3^$s$s3" . $rem;
          $next =~ s/\/[^\/]+$/$diff/;
          $next =~ s/($cand)([@positions]_\d\^\d.*)($cand)[@positions]_\d\^\d/$1$2/;
          push(@tmp,$next);
        } elsif ( $cnt =~ /$cand([@positions])_(\d)\^(\d)/ ) {
          my $y = $1;
          my $c2 = $2;
          my $s2 = $3;
          my $next = $cnt;

          my $test = $next;
          my $sign = 1;
          $test =~ s/$cand(.+)$//;
          while ( $test =~ /D/ ) {
            $test =~ s/D//;
            $sign *= -1;
          }
          $next = " -" . $next if ( $sign == -1 );

          my $diff = "S$q($x-$y)_$c$c2^$s$s2" . $rem;
          $next =~ s/\/[^\/]+$/$diff/;
          $next =~ s/$cand([@positions])_\d\^\d//;

          push(@tmp,$next);
        }
      }
    }# $cnt
    @cont = @tmp;
    &show_contract([@cont]) unless ( $show_flag == 0 );
    $ctr++;
    last if ( $ctr == $this->{num_f} );
    $flag = 0;
    foreach my $con (@cont) {
      $flag = 1 if ( $con =~ /\// );
    }
  }# while

  @cont = &signManage(@cont);
#  &show_contract([@cont],"# Contraction of$this->{corr}:");
  @cont = $this->reorder($type_label,@cont);
  &show_contract([@cont],"# Contraction of$this->{corr}:");
#       unless ( $show_flag == 0 );
  $this->{tracef} = [@cont];
  return @cont;
}

sub reorder {
  my $this = shift;
  my $type_label = shift;
  my @cont = @_;
  my @result = ();
  my @gamma = @{$this->{gamma}};
  foreach my $term (@cont) {
    $term =~ s/^\s+//;
    $term =~ s/\s+$//;
    my $newt = " ";
    my @tmp = split('\s+',$term);
    my $sign = shift(@tmp);
    my $coef = $this->{coef}->[0];
    if ( $sign eq "-" ) {
      $coef *= -1;
    }
    if ( $coef == 1 ) {
      $coef = "+";
    } elsif ( $coef == -1 ) {
      $coef = "-";
    } elsif ( $coef > 0 ) {
      $coef = "+$coef";
    }
    $newt .= $coef . " ";
    my $col = 1;
    my $flag = 2;
    my $tmpc = " ";
    while( $#tmp > -1 ) {
      my @tmp2 = ();
      foreach my $prop (@tmp) {
        if ( $prop =~ /_$col(\d)\^(\d)/ ) {
          $col = $1;
          my $spin = $2;
          foreach my $g (@gamma) {
            if ( $g =~ /\^$spin$/ ) {
              $tmpc .= "$g ";
              last;
            }
          }
          $tmpc .= "$prop ";
        } else {
          push(@tmp2,$prop);
        }
      }
      if ( $#tmp == $#tmp2 || $#tmp2 == -1 ) {
        $col = $flag;
        $flag++;
        $newt .= "Tr[$tmpc] " unless ( $tmpc eq " " );
        $tmpc = " ";
      }
      @tmp = @tmp2;
    }
    $newt =~ s/\^\d+//g;
    $newt =~ s/_\d+//g;
    if ( $type_label < 2 ) {
      next if ( $newt =~ /x-x/ && $type_label == 1 );# neglect qloop diagrams
      next unless ( $newt =~ /x-x/ || $type_label == 1 );# neglect qloop diagrams
    }
    push(@result,$newt);
  }
  return @result;
}

sub show_contract {
  my $cont_r = shift;
  my $comment = ( shift || "" );
  print "$comment\n" unless ( $comment eq "" );
  foreach my $con (@{$cont_r}) {
    print "$con\n";
  }
  print "\n";
  return;
}

sub show_object {
  my $this = shift;
  my $comment = ( shift || "" );
  print "$comment\n" unless ( $comment eq "" );
#  print "$this->{corr}\n";
  my @coef = @{$this->{coef}};
  die unless ( $coef[1] == 0 );# correlators supposed to be real
  print sprintf("Coef: Re = %f, Im = %f\n",@coef);
  my @gamma = @{$this->{gamma}};
#  print "Spin Structure: @gamma\n\n";
  return;
}

sub Hconj {
  # Create a hermitian conjugate of an operator
  #
  # * this hermitian conjugate is in Minkowski
  # but Euclidean correlator is supposed to be calculated
  # with this hermitian conjugate
  #
  my $O = shift;
  $O =~ s/^\s+//;
  my @tmp = split('\s+',$O);
  my $nic = 1;
  my $nis = 1;
  my $niv = 1;
  foreach my $part (@tmp) {
    if ( $part =~ /_(\d+)/ ) {# color index
      my @ids = split('\s*',$1);
      $nic = &Functions::max($nic-1,@ids) + 1;
    }
    if ( $part =~ /\^(\d+)/ ) {
      my @ids = split('\s*',$1);
      $nis = &Functions::max($nis-1,@ids) + 1;
    }
    if ( $part =~ /mu(\d+)/ ) {
      my @ids = split('\s*',$1);
      $niv = &Functions::max($niv-1,@ids) + 1;
    }
  }

  my $Hc = "";
  my $sign = 1;
  my @spins = ();
  my @colors = ();
  my @vectors = ();
  my @indicator_c = ();
  my @indicator_s = ();
  my @indicator_v = ();
  until ( $#tmp == -1 ) {
    my $part = pop(@tmp);
    if ( $part =~ /([@quarks])v_(\d+)\^(\d+)/ ) {
      my $ic = &Functions::inverse_array($2,[@indicator_c]);
      if ( $ic == -1 ) {
        $Hc .= " $1" . "b_$nic";
        push(@colors,$nic);
        push(@indicator_c,$2);
        $nic++;
      } else {
        $Hc .= " $1" . "b_$colors[$ic]";
      }

      my $is = &Functions::inverse_array($3,[@indicator_s]);
      if ( $is == -1 ) {
        $Hc .= "^$nis";
        push(@spins,$nis);
        push(@indicator_s,$3);
        $nis++;
      } else {
        $Hc .= "^$spins[$ic]";
      }
    } elsif ( $part =~ /([@quarks])b_(\d+)\^(\d+)/ ) {
      my $ic = &Functions::inverse_array($2,[@indicator_c]);
      if ( $ic == -1 ) {
        $Hc .= " $1" . "v_$nic";
        push(@colors,$nic);
        push(@indicator_c,$2);
        $nic++;
      } else {
        $Hc .= " $1" . "v_$colors[$ic]";
      }

      my $is = &Functions::inverse_array($3,[@indicator_s]);
      if ( $is == -1 ) {
        $Hc .= "^$nis";
        push(@spins,$nis);
        push(@indicator_s,$3);
        $nis++;
      } else {
        $Hc .= "^$spins[$is]";
      }
    } elsif ( $part eq "i" ) {
      $Hc .= " - i";
    } elsif ( $part =~ /^\s*g([^\^]+)\^(\d+)\^(\d+)/ ) {
      my $gam = $1;
      my $s1 = $2;
      my $s2 = $3;
      if ( $gam eq "5" ) {
        $Hc .= " - g$1";
      } elsif ( $1 =~ /mu(\d+)$/ ) {
        my $imu = &Functions::inverse_array_char($1,[@indicator_v]);
        if ( $imu == -1 ) {
          $Hc .= " - gmu$niv";
          push(@vectors,$niv);
          push(@indicator_v,$1);
          $niv++;
        } else {
          $Hc .= " - gmu$vectors[$imu]";
        }
      } elsif ( $gam =~ /^mu(\d+)([L,R])$/ ) {
        my $imu = &Functions::inverse_array_char($1,[@indicator_v]);
        if ( $imu == -1 ) {
          $Hc .= " gmu$niv$2";
          push(@vectors,$niv);
          push(@indicator_v,$1);
          $niv++;
        } else {
          $Hc .= " gmu$vectors[$imu]$2";
        }
      }

      my $is = &Functions::inverse_array($s2,[@indicator_s]);
      if ( $is == -1 ) {
        $Hc .= "^$nis";
        push(@spins,$nis);
        push(@indicator_s,$s2);
        $nis++;
      } else {
        $Hc .= "^$spins[$is]";
      }

      $is = &Functions::inverse_array($s1,[@indicator_s]);
      if ( $is == -1 ) {
        $Hc .= "^$nis";
        push(@spins,$nis);
        push(@indicator_s,$s1);
        $nis++;
      } else {
        $Hc .= "^$spins[$is]";
      }
    }
  }
  return $Hc;
}

1;
