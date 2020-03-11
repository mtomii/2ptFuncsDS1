package Functions;
use strict;
use warnings;
#use Math::Round;
use File::Path;
use File::Path;
use Math::Trig 'pi';
use Math::Trig 'acos';

our $pi = pi;

sub get_pi {return $pi}

sub cube_root {
  my $val = shift;
  return $val / abs($val) * abs($val)**(1.0/3.0);
}

sub min {# get minimum value in input array
  my $min = shift;
  foreach(@_) {
    $min = $_ if ( $min > $_ );
  }

  return $min;
}

sub max {# get max value in input array
  my $max = shift;
  foreach(@_) {
    $max = $_ if ( $max < $_ );
  }
  return $max;
}

sub make_path {
  my $target_dir = shift;
  eval{
    mkpath($target_dir, {verbose => 1}, 0755);
  };
  if ( $@ ) {
    die "Failed to make $target_dir\n";
  }
  return;
}

sub randn {
  my ($m, $sigma) = @_;
  my ($r1, $r2) = (rand(), rand());
  while ($r1 == 0) { $r1 = rand(); }
  return ($sigma * sqrt(-2 * log($r1)) * sin(2 * pi * $r2)) + $m;
}

sub get_time_stump {
  my $file = shift;
  my $modtime = 0;
  $modtime = (stat($file))[9] if (-e $file);
  return $modtime;
}

sub time_stump {
  my $file = shift;
  my $modtime = &get_time_stump($file);
  my $sec;
  my $min;
  my $hour;
  my $mday;
  my $mon;
  my $year;
  ($sec, $min, $hour, $mday, $mon, $year) = localtime($modtime);
  $year += 1900;
  $mon++;
  print "time stump<$file>  $year.$mon.$mday $hour:" .
      sprintf("%02d:%02d\n",$min,$sec);
  return;
}

sub get_unexist_or_oldest_file {
  my @files = @_;
  my @stump = (59,59,23,31,11,999);
  my @result = ();
  my $Wflag = 0;
  foreach my $file(@files) {
    my $modtime = &get_time_stump($file);
    my @stump1 = localtime($modtime);
    for ( my $id = $#stump ; $id >= 0 ; $id-- ) {
      last if ( $stump1[$id] > $stump[$id]);
      if ($stump1[$id] < $stump[$id]) {
	@stump = @stump1;
	@result = ($file);
	last;
      }
      push(@result,$file) if ( $id == 0 );
    }
  }
  return @result;
}

sub get_newest_file {
  my @files = @_;
  my @stump = (0,0,0,0,0,0);
  my @result = ();
  my $Wflag = 0;
  foreach my $file(@files) {
    my $modtime = &get_time_stump($file);
    my @stump1 = localtime($modtime);
    for ( my $id = $#stump ; $id >= 0 ; $id-- ) {
      last if ( $stump1[$id] < $stump[$id]);
      if ($stump1[$id] > $stump[$id]) {
	@stump = @stump1;
	@result = ($file);
	last;
      }
      push(@result,$file) if ( $id == 0 );
    }
  }
  return @result;
}

sub gauss_bracket {
  my $val = shift;
  return int($val) - 1 if ( $val <= 0.0 );
  return int($val);
}

sub factorial {
  my $n = shift;
  my $f = 1;
  return $f if ( $n <= 1 );
  for ( my $i = 2 ; $i <= $n ; $i++ ) {
    $f *= $i;
  }
  return $f;
}

sub log10 {
  my $x = shift;
  return log($x) / log(10);
}

sub cos_deg {
  my $deg = shift;
  my $rad = $deg * pi / 180.0;
  return cos($rad);
}

sub zeta_func_old {
  # zeta(n) = \sum_{k=1}^\infty k^{-n}
  # input1 : n
  # input2 : <option> : stoping condition
  my $int = shift;
  my $precision = ( shift || 1.0E-16 );
  my $flag = 1;
  my $zeta = 0;
  my $k = 1;
  return $pi**2 / 6 if ( $int == 2 );
  return $pi**4 / 90 if ( $int == 4 );
  return $pi**6 / 945 if ( $int == 6 );
  return "infinity" if ( $int <  2 );
  $precision = 1.0E-14 if ( $int < 2.3 && $precision < 1.0E-14 );
  while($flag==1) {
    $zeta += 1/$k**$int;
#p    print "tmp: $zeta\n";
    $flag = 0 if ( $k**$int * $precision > 1 );
    $k++;
  }
  return $zeta;
}

sub new_int_array {
  my $num = shift;
  my $init = ( shift || 0 );
  my @array = ();
  for( my $i = 0 ; $i <= $num ; $i++ ) {
    push(@array, $init);
  }

  return @array;
}

sub new_float_array {
  my $num = shift;
  my $init = ( shift || 0.0 );
  my @array = ();
  for( my $i = 0 ; $i <= $num ; $i++ ) {
    push(@array, $init);
  }
  return @array;
}

sub new_char_array {
  my $num = shift;
  my $init = ( shift || "" );
  my @array = ();
  for( my $i = 0 ; $i < $num ; $i++ ) {
    push(@array, $init);
  }

  return @array;
}

sub new_ref_array {
  my $num = shift;
  my $init = ( shift || [] );
  my @array = ();
  for( my $i = 0 ; $i <= $num ; $i++ ) {
    push(@array, $init);
  }
  return @array;
}

sub new_int_mtx {
  my $num1 = shift;
  my $num2 = ( shift || $num1 );
  my $init = ( shift || 0 );
  my @mtx = ();

  for ( my $i = 0 ; $i <= $num1 ; $i++ ) {
    my @array = &new_int_array($num2, $init);
    push @mtx, [@array];
  }

  return [@mtx];
}

sub new_float_mtx {
  my $num1 = shift;
  my $num2 = shift;
  my $init = ( shift || 0.0 );
  my @mtx = ();
  for ( my $i = 0 ; $i <= $num1 ; $i++ ) {
    my @array = &new_float_array($num2, $init);
    push @mtx, [@array];
  }

  return [@mtx];
}

sub get_id_mtx {
  my $size = shift;
  my @id = @{&new_float_mtx($size,$size)};
  for ( my $i = 0 ; $i <= $size ; $i++ ) {
    $id[$i]->[$i] = 1.0;
  }
  return @id;
}

sub new_char_mtx {
  my $num1 = shift;
  my $num2 = ( shift || $num1 );
  my $init = ( shift || "" );
  my @mtx = ();

  for ( my $i = 0 ; $i <= $num1 ; $i++ ) {
    my @array = &new_char_array($num2, $init);
    push @mtx, [@array];
  }

  return [@mtx];
}

sub new_ref_mtx {
  my $num1 = shift;
  my $num2 = ( shift || $num1 );
  my $init = ( shift || [] );
  my @mtx = ();

  for ( my $i = 0 ; $i <= $num1 ; $i++ ) {
    my @array = &new_ref_array($num2, $init);
    push @mtx, [@array];
  }

  return [@mtx];
}

sub check_existence_numeric_array {
  # check existence in input numeric array
  my $num = shift;
  my $ref = shift;
  my @array = @$ref;
  my $existence = 0;

  foreach my $element (@array) {
    $existence = 1 if ( $element == $num );
  }

  return $existence;
}

sub add_array {
  my ( $ary, $array ) = @_;
  my @a1 = @$ary;
  my @a2 = @$array;

  if ( $#a1 < $#a2 ) {
    @a1 = @$array;
    @a2 = @$ary;
  }

  print "Warning : subroutine add_array : the size of two arrays are not equal to each other( $#a1 : $#a2 )\n" if ( $#a1 != $#a2 );

  my @result = @a1;

  for ( my $i = 0 ; $i <= $#a2 ; $i++ ) {
    $result[$i] += $a2[$i];
  }

  return @result;
}

sub add_mtx {
  my $m1_ref = shift;
  my $m2_ref = shift;

  my @m1 = @$m1_ref;
  my @m2 = @$m2_ref;
  {
    my @size1 = &get_size_mtx(@m1);
    my @size2 = &get_size_mtx(@m2);
    if ( $size1[0] == -1 || $size2[0] == -1 ) {
      return @m1 unless ( $size1[0] == -1 );
      return @m2 unless ( $size2[0] == -1 );
      return ();
    } elsif ( $size1[1] == -1 || $size2[1] == -1 ) {
      return @m1 unless ( $size1[1] == -1 );
      return @m2;
    } else {
      die "dead at sub: add_mtx A @size1 : @size2\n"
	  unless ( $size1[0] == $size2[0] && $size1[1] == $size2[1] );
    }
  }

  my @result = ();
  for ( my $i = 0 ; $i <= $#m1 ; $i++ ) {
    my @m1i = @{$m1[$i]};
    my @m2i = @{$m2[$i]};
    return () if ( $#m1i == -1 && $#m2i == -1 );
    return (&divide_mtx([@m2],-1.0)) if ( $#m1i == -1 );
    return @m1 if ( $#m2i == -1 );
    die "dead at sub: add_mtx B,$i, $#m1i $#m2i\n" unless ($#m1i == $#m2i);
    my @mi = &add_array([@m1i],[@m2i]);
    push @result, [@mi];
  }
  return @result;;
}

sub sub_array {
  my ( $ary, $array ) = @_;
  my @a1 = @$ary;
  my @a2 = @$array;

  my $n = $#a2;
  $n = $#a1 if ( $n > $#a1 );
  print "Warning : subroutine sub_array : the size of two arrays are not equal to each other( $#a1 : $#a2 )\n" if ( $#a1 != $#a2 );

  my @result = @a1;
  for ( my $i = 0 ; $i <= $#a2 ; $i++ ) {
    $result[$i] -= $a2[$i];
  }

  return @result;
}

sub sub_mtx {
  my $m1_ref = (shift||[]);
  my $m2_ref = (shift||[]);

  my @m1 = @$m1_ref;
  my @m2 = @$m2_ref;
  {
    my @size1 = &get_size_mtx(@m1);
    my @size2 = &get_size_mtx(@m2);
    if ( $size1[0] == -1 || $size2[0] == -1 ) {
      return @m1 unless ( $size1[0] == -1 );
      return @m2 unless ( $size2[0] == -1 );
      return ();
    } elsif ( $size1[1] == -1 || $size2[1] == -1 ) {
      return @m1 unless ( $size1[1] == -1 );
      return &divide_mtx([@m2],-1.0);
    } else {
      die "dead at sub: sub_mtx A @size1 : @size2\n"
	  unless ( $size1[0] == $size2[0] && $size1[1] == $size2[1] );
    }
  }

  my @result = ();
  for ( my $i = 0 ; $i <= $#m1 ; $i++ ) {
    my @m1i = @{$m1[$i]};
    my @m2i = @{$m2[$i]};
    return () if ( $#m1i == -1 && $#m2i == -1 );
    return (&divide_mtx([@m2],-1.0)) if ( $#m1i == -1 );
    return @m1 if ( $#m2i == -1 );
    die "dead at sub: sub_mtx B,$i, $#m1i $#m2i\n" unless ($#m1i == $#m2i);
    my @mi = &sub_array([@m1i],[@m2i]);
    push @result, [@mi];
  }

  return @result;;
}

sub add_array_sq {
  my ( $ary, $array ) = @_;
  my @a1 = @$ary;
  my @a2 = @$array;

  if ( $#a1 != $#a2 ) {
    print "Error : subroutine sum_array_sq : the size of two arrays are not equal to each other.\n";
    die;
  }

  my @result = @a1;

  for ( my $i = 0 ; $i <= $#a2 ; $i++ ) {
    $result[$i] += $a2[$i] **2;
  }

  return @result;
}

sub add_mtx_sq {
  my $m1_ref = shift;
  my $m2_ref = shift;
  my @m1 = @$m1_ref;
  my @m2 = @$m2_ref;
  die "dead at sub: add_mtx_sq A $#m1 $#m2\n" unless ($#m1 == $#m2);

  my @result = ();
  for ( my $i = 0 ; $i <= $#m1 ; $i++ ) {
    my @m1i = @{$m1[$i]};
    my @m2i = @{$m2[$i]};
    die "dead at sub: add_mtx_sq B,$i, $#m1i $#m2i\n" unless ($#m1i == $#m2i);
    my @mi = &add_array_sq([@m1i],[@m2i]);
    push @result, [@mi];
  }

  return @result;
}

sub sub_array_sq {
  my ( $ary, $array ) = @_;
  my @a1 = @$ary;
  my @a2 = @$array;

  if ( $#a1 != $#a2 ) {
    print "Error : subroutine sub_array_sq : the size of two arrays are not equal to each other.\n";
    die;
  }

  my @result = @a1;

  for ( my $i = 0 ; $i <= $#a2 ; $i++ ) {
    $result[$i] -= $a2[$i] **2;
  }

  return @result;
}

sub divide_array {
  my $ref = shift;
  my $divider = shift;
  my @array = @$ref;

  for ( my $i = 0 ; $i <= $#array ; $i++ ) {
    $array[$i] /= $divider;
  }

  return @array;
}

sub pow_array_comp {
  my $ref = shift;
  my $alp = shift;

  my @result = ();
  foreach my $mi (@$ref) {
    push (@result, $mi**$alp);
  }

  return @result;
}

sub divide_mtx {
  my $ref = shift;
  my $divider = shift;

  my @result = ();
  foreach my $mi (@$ref) {
    my @tmp = ();
    foreach my $mij (@$mi) {
      push (@tmp,$mij/$divider);
    }
    push @result, [@tmp];
  }

  return @result;
}

sub divide_mtx_each_comp {
  my $ref = shift;
  my $divider_ref = shift;

  my @result = ();
  my @divider = @$divider_ref;
  for ( my $i = 0 ; $i <= $#divider ; $i++ ) {
    my @tmp = ();
    my @nume = @{$ref->[$i]};
    for ( my $j = 0 ; $j <= $#nume ; $j++ ) {
      push(@tmp,$nume[$j]/$divider[$i]->[$j]);
    }
    push @result, [@tmp];
  }

  return @result;
}

sub mult_mtx {
  my $ref = shift;
  my $multer = (shift||0);

  my @result = ();
  foreach my $mi (@$ref) {
    my @tmp = ();
    foreach my $mij (@$mi) {
      push (@tmp,$mij*$multer);
    }
    push @result, [@tmp];
  }

  return @result;
}

sub pow_mtx_comp {
  my $ref = shift;
  my $alp = shift;

  my @result = ();
  foreach my $mi (@$ref) {
    my @tmp = ();
    foreach my $mij (@$mi) {
      push (@tmp,$mij**$alp);
    }
    push @result, [@tmp];
  }

  return @result;
}

sub transpose_mtx {
  my @m = @_;

  my $n1 = $#m;
  return () if ( $n1 == -1 );
  return () unless ( defined $m[0] );
  my $n2;
  {
    my @m0 = @{$m[0]};
    $n2 = $#m0;
  }

  my @result = ();
  for ( my $i = 0 ; $i <= $n2 ; $i++ ) {
    my @tmp = ();
    for ( my $j = 0 ; $j <= $n1 ; $j++ ) {
      push(@tmp,$m[$j]->[$i]);
    }
    push @result, [@tmp];
  }

  return @result;
}

sub get_size_mtx {
  my @mtx = @_;
  my $m1 = $#mtx;
  return (-1,-1) if ( $m1 == -1 );
  my @mtx0 = @{$mtx[0]};
  my $m2 = $#mtx0;
  return ( $m1, $m2 );
}

sub mtx_mtx {
  my $m1_ref = ( shift || [] );
  my $m2_ref = ( shift || [] );

  my @m1 = @$m1_ref;
  return () if ( $#m1 == -1 );
  my @size1 = &get_size_mtx(@m1);
  my @m2 = @$m2_ref;
  @m2 = @{&new_float_mtx($size1[1],0,0)} if ( $#m2 == -1 );
  my @size2 = &get_size_mtx(@m2);
#  print "@size1 : @size2\n";
  die "sizes of input matrices are not mutch! @size1 : @size2\n"
      unless ( $size1[1] == $size2[0] );

  my $s1 = $size1[0];
  my $s2 = $size2[1];
  my $s = $size1[1];
  my @result = @{&new_float_mtx($s1,$s2)};
  for ( my $i = 0 ; $i <= $s1 ; $i++ ) {
    for ( my $j = 0 ; $j <= $s2 ; $j++ ) {
      for ( my $k = 0 ; $k <= $s ; $k++ ) {
	$result[$i]->[$j] += $m1[$i]->[$k] * $m2[$k]->[$j];
      }
    }
  }
  return @result;
}

sub basic_trf_exch {
  my $mtx_ref = shift;
  my $i = shift;
  my $j = shift;
  my @mtx = @$mtx_ref;

  my $tmp_ref = $mtx[$i];
  $mtx[$i] = $mtx[$j];
  $mtx[$j] = $tmp_ref;

  return @mtx;
}

sub basic_trf_norm {
  my $mtx_ref = shift;
  my $i = shift;
  my $norm = ( shift || $mtx_ref->[$i]->[$i] );
  my @mtx = @$mtx_ref;

  $mtx[$i] = [&divide_array($mtx[$i],$norm)];
  return @mtx;
}

sub basic_trf_comb {
#      @mtx = &basic_trf_comb([@mtx],$n,$i,$ratio);
  my $mtx_ref = shift;
  my $i = shift;
  my $j = shift;
  my $ratio = ( shift || $mtx_ref->[$j]->[$i] );
  my @mtx = @$mtx_ref;

  my @tmp = &mult_array($mtx[$i],$ratio);
  $mtx[$j] = [&sub_array($mtx[$j],[@tmp])];

  return @mtx;
}

sub inv_mtx {
  my $mtx_ref = ( shift || [] );
  my $trg_ref = ( shift || [] );
  my @mtx = @$mtx_ref;
  return () if ( $#mtx == -1 );
  my @size_mtx = &get_size_mtx(@mtx);
  die "Not squared matrix: $size_mtx[0] x $size_mtx[1]\n"
      unless ( $size_mtx[0] == $size_mtx[1] );
  my $size = $size_mtx[0];

  my @trg = @$trg_ref;
  return () if ( $#trg == -1 );
  for ( my $n = 0 ; $n <= $size ; $n++ ) {
    if ( $mtx[$n]->[$n] == 0.0 ) {
      my $i = 0;
      for ( my $k = $n+1 ; $k <= $size ; $k++ ) {
	unless ( $mtx[$k]->[$k] == 0.0 ) {
	  $i = $k;
	  last;
	}
      }
      die "Singular matrix!\n" if ( $i == 0 );
      @mtx = &basic_trf_exch([@mtx],$n,$i);
      @trg = &basic_trf_exch([@trg],$n,$i);
    }

    my $norm = $mtx[$n]->[$n];
    @mtx = &basic_trf_norm([@mtx],$n,$norm);
    @trg = &basic_trf_norm([@trg],$n,$norm);

    for ( my $i = 0 ; $i <= $size ; $i++ ) {
      next if ( $i == $n);
      my $ratio = $mtx[$i]->[$n];
      unless ( $ratio == 0.0 ) {
	@mtx = &basic_trf_comb([@mtx],$n,$i,$ratio);
	@trg = &basic_trf_comb([@trg],$n,$i,$ratio);
      }
    }
  }# $n
  return @trg;
}

sub inv_mtx_detailed {
  my $mtx_ref = ( shift || [] );
  my $trg_ref = ( shift || [] );
  my @mtx = @$mtx_ref;
  return () if ( $#mtx == -1 );
  my @size_mtx = &get_size_mtx(@mtx);
  die "Not squared matrix: $size_mtx[0] x $size_mtx[1]\n"
      unless ( $size_mtx[0] == $size_mtx[1] );
  my $size = $size_mtx[0];

  my @trg = @$trg_ref;
  return () if ( $#trg == -1 );
  print "Start\n";
  &show_mtx(@mtx);
  print "/\n";
  &show_mtx(@trg);
  print "\n";
  for ( my $n = 0 ; $n <= $size ; $n++ ) {
    if ( $mtx[$n]->[$n] == 0.0 ) {
      my $i = 0;
      for ( my $k = $n+1 ; $k <= $size ; $k++ ) {
	unless ( $mtx[$k]->[$k] == 0.0 ) {
	  $i = $k;
	  last;
	}
      }
      die "Singular matrix!\n" if ( $i == 0 );
      @mtx = &basic_trf_exch([@mtx],$n,$i);
      @trg = &basic_trf_exch([@trg],$n,$i);
    }

    my $norm = $mtx[$n]->[$n];
    @mtx = &basic_trf_norm([@mtx],$n,$norm);
    @trg = &basic_trf_norm([@trg],$n,$norm);
    print "Normalized\n";
    &show_mtx(@mtx);
  print "/\n";
    &show_mtx(@trg);
    print "\n";

    for ( my $i = 0 ; $i <= $size ; $i++ ) {
      next if ( $i == $n);
      my $ratio = $mtx[$i]->[$n];
      unless ( $ratio == 0.0 ) {
	@mtx = &basic_trf_comb([@mtx],$n,$i,$ratio);
	@trg = &basic_trf_comb([@trg],$n,$i,$ratio);
      }
    }
    print "Swiped\n";
    &show_mtx(@mtx);
  print "/\n";
    &show_mtx(@trg);
    print "\n";
  }# $n
  &show_mtx(@mtx);
  return @trg;
}

sub inv_mtx_id {
  my $mtx_ref = ( shift || [] );
  my @mtx = @$mtx_ref;
  return &inv_mtx($mtx_ref,[&get_id_mtx($#mtx)]);
}

sub inv_mtx_id_detailed {
  my $mtx_ref = ( shift || [] );
  my @mtx = @$mtx_ref;
  return &inv_mtx_detailed($mtx_ref,[&get_id_mtx($#mtx)]);
}

sub show_mtx {
  my @m = @_;
  foreach my $mi_ref (@m) {
    my @mi = @$mi_ref;
    print "@mi\n";
  }
  return;
}

sub write_mtx {
  my $m = shift;
  my $handle = shift;
  foreach my $mi_ref (@$m) {
    my @mi = @$mi_ref;
    print $handle "@mi\n";
  }
  return;
}

sub complex_pow {
  my $in = shift;
  my $pow = shift;

  my $x = $in->[0];
  my $y = $in->[1];
  my $r = sqrt($x*$x + $y*$y);

  return (0,0) if ($r==0);
  my $theta = acos($x/$r);
  $theta *= -1 if ( $y < 0 );
#  $theta = 2*pi-$theta if ( $y < 0 );

  $r = $r**$pow;
  $theta *= $pow;
  return ($r*cos($theta),$r*sin($theta));
}

sub complex_sum {
  my $in1 = shift;
  my $in2 = shift;
  return ($in1->[0]+$in2->[0],$in1->[1]+$in2->[1]);
}

sub complex_comb {
  my $in1 = shift;
  my $in2 = shift;
  my $c1 = shift;
  my $c2 = shift;
  return ($c1*$in1->[0]+$c2*$in2->[0],$c1*$in1->[1]+$c2*$in2->[1]);
}

sub complex_mult {
  my $in1 = shift;
  my $in2 = shift;
  return ($in1->[0]*$in2->[0] - $in1->[1]*$in2->[1],
          $in1->[0]*$in2->[1] + $in1->[1]*$in2->[0]);
}

sub complex_div {
  my $in1 = shift;
  my $in2 = shift;

  my $r2_sq = $in2->[0]**2 + $in2->[1]**2;
  return (($in1->[0]*$in2->[0]+$in1->[1]*$in2->[1])/$r2_sq,
         -($in1->[0]*$in2->[1]-$in1->[1]*$in2->[0])/$r2_sq);
}

sub complex_abs {
  my @in = @_;
  return sqrt( $in[0]*$in[0] + $in[1]*$in[1] );
}

sub mult_array {
  my $ref = shift;
  my $multer = shift;
  my @array = @$ref;

  for ( my $i = 0 ; $i <= $#array ; $i++ ) {
    $array[$i] *= $multer;
  }
  return @array;
}

sub divide_array_each_comp {
  my $ref = shift;
  my $divider_ref = shift;
  my @divider = @$divider_ref;

  my @result = ();
  for ( my $i = 0 ; $i <= $#divider ; $i++ ) {
    push (@result, $ref->[$i]/$divider[$i]);
  }
  return @result;
}

sub mult_array_each_comp {
  my $ref = shift;
  my $multer_ref = shift;
  my @multer = @$multer_ref;

  my @result = ();
  for ( my $i = 0 ; $i <= $#multer ; $i++ ) {
    push (@result, $ref->[$i]*$multer[$i]);
  }
  return @result;
}

sub sqrt_array {
  my @in = @_;
  my @result = ();

  foreach my $comp (@in) {
    push(@result,sqrt($comp));
  }

  @result;
}

sub sqrt_or0_array {
  my @in = @_;
  my @result = ();
  foreach my $comp (@in) {
    push(@result,sqrt($comp)) if ( $comp > 0 );
    push(@result,0) if ( $comp <= 0 );
  }
  @result;
}

sub scalar_prod {
  my $v1_ref = shift;
  my $v2_ref = shift;

  my @v1 = @$v1_ref;
  my @v2 = @$v2_ref;

  my $prod = 0.0;
  my $n = $#v2;
  $n = $#v1 if ( $n > $#v1 );
  print "Warning : subroutine scalar_prod : the size of two arrays are not equal to each other( $#v1 : $#v2 )\n" if ( $#v1 != $#v2 );

  for ( my $i = 0 ; $i <= $n ; $i++ ) {
    $prod += $v1[$i] * $v2[$i];
  }

  return $prod;
}

sub metric_prod {
  my $v1_ref = shift;
  my $mtx_ref = shift;
  my $v2_ref = shift;

  my @v1 = @$v1_ref;
  my @v2 = @$v2_ref;

  my $prod = 0.0;

  for ( my $i = 0 ; $i <= $#v1 ; $i++ ) {
    for ( my $j = 0 ; $j <= $#v2 ; $j++ ) {
      $prod += $v1[$i] * $mtx_ref->[$i]->[$j] * $v2[$j];
    }
  }

  return $prod;
}

sub norm_array {
  my $v_ref = shift;
  return sqrt(&scalar_prod($v_ref,$v_ref));
}

sub add_count_array {
  my ( $d, $c ) = @_;
  my @dat = @$d;
  my @cnt = @$c;

  my $num = $#dat;
  if ( $#dat != $#cnt ) {
    print "Warning : subroutine add_count_array : the size of two arrays are not equal to each other.\n";
    $num = $#cnt if ( $#dat > $#cnt );
  }

  for ( my $i = 0 ; $i <= $num ; $i++ ) {
    $cnt[$i]++ if ( $dat[$i] != 0.0 );
  }

  return @cnt;
}

sub get_combination {
  # Get nCm
  my $n = shift;
  my $m = shift;
  return &factorial($n)/&factorial($m)/&factorial($n-$m);
}

sub get_permutation {
  # input  : numeric array
  # output : number of possible permutations
  my @array = @_;

  my @tmp = ();
  my @nums = ();
  foreach my $val (@array) {
    my $idx = &inverse_array($val,[@nums]);
    if ( $idx == -1 ) {
      push(@tmp,1);
      push(@nums,$val);
    } else {
      $tmp[$idx]++;
    }
  }
  my $p = 1;
  foreach my $n (@tmp) {
    $p *= &factorial($n);
  }
  my $f = &factorial( $#array + 1 );

  return $f / $p;
}

sub get_all_combination_array {
# if you call this with input number and numeric array,
#        you can get all combination
  my $num   = shift;
  my @array = sort {$a <=> $b} @_;

  my @ary = ($array[0]);
  my $j = 0;
  for ( my $i = 1 ; $i <= $#array ; $i++ ) {
    next if ( $array[$i] == $ary[$j] );
    push ( @ary, $array[$i] );
    $j++;
  }

  print "Functions::get_combination : the input number($num) is larger than the size of sorted input array($j)\n" if ( $num >= $j );

  return @array;
}

sub reorder_ascend {
  my @array = sort {$a <=> $b} @_;
  return @array;
}

sub inverse_array {
  # If you call this subroutine with input value and ref(array): ($val,$ref),
  # you can get the index $idx, which is such that
  # @$ref[$idx] = $val
  my $val = shift;
  my $ref = shift;
  my $idx = -1;
  my $tmp = 0;

  foreach my $component (@$ref) {
    if( $component == $val ) {
      $idx = $tmp;
      last;
    }
    $tmp++;
  }
  return $idx;
}

sub inverse_array_rough {
  # If you call this subroutine with input value and ref(array): ($val,$ref),
  # you can get the index $idx, which is such that
  # @$ref[$idx] = $val
  my $val = shift;
  my $ref = shift;
  my $idx = -1;
  my $tmp = 0;

  foreach my $component (@$ref) {
    if( $component+$val == 0 ) {
      if ( $component == 0 ) {
        $idx =$tmp;
        last;
      }
      next;
    }
#    print abs($component-$val) . " / " . abs($component+$val) . " = " . (abs($component-$val)/abs($component+$val)) . "\n";
    if( abs($component-$val)/abs($component+$val) < 1e-06 ) {
      $idx = $tmp;
      last;
    }
    $tmp++;
  }
  return $idx;
}

sub inverse_array_char {
  # charactor array version of "sub inverse_array"
  my $char = shift;
  my $ref  = shift;
  my $idx  = -1;
  my $tmp  = 0;

#p  print "SSS$char\n";
  foreach my $component (@$ref) {
#p    print "$char\n";
    if ( $component eq $char ) {
      $idx = $tmp;
#p      print "Hit!\n" . "$component\n";
      last;
    }
#p    print "$component\n\n\n";
    $tmp++;
  }
  return $idx;
}

sub cmplx_sqrt {
  my @z = @_;
  my $abs = &norm_array([@z]);
  my $cos_arg_z = $z[0]/$abs;

  my @phase_part = (sqrt((1.0+$cos_arg_z)/2.0),
		    sqrt((1.0-$cos_arg_z)/2.0));
  $phase_part[0] *= -1.0 if ( $z[1] < 0.0 );

  my @sqrt_z = &mult_array([@phase_part],sqrt($abs));
  return @sqrt_z;
}

sub cmplx_log {
  my @z = @_;
  my $abs = &norm_array([@z]);
  my $cos_arg_z = $z[0]/$abs;
  die "sub:cmplx_log:: argument is 0, cannot take log\n" if ($abs==0.0);
  my $re = log($abs);
  my $im = acos($cos_arg_z);

  return ($re,$im);
}

sub get_jks_single {
  # input1: data_ref: [$conf1,$conf2,...,$confN]
  # input2: binsize
  # if you call this subroutine, you can get jackknif sample
  # as the form ($Conf1,$Conf2,...,$Conf{'N/binsize'})
  my $data_ref = shift;
  my $bin_size = (shift||1);
  my @data = @$data_ref;

  my $total = 0;
  foreach my $cdata ( @data ) {
    $total += $cdata;
  }

  my @jks = ();
  my $flag = 0;
  my $cnt = 0;
  while ( $flag == 0 ) {
    my $tmp = $total;
    my $ctr = 0;
    for ( my $ic = $cnt * $bin_size ; $ic < ($cnt+1) * $bin_size ; $ic++ ) {
      $tmp -= $data[$ic];
      $ctr++;
      if ( $ic == $#data ) {
	$flag = 1;
	last;
      }
    }
    $tmp /= $#data+1-$ctr;
    push (@jks,$tmp);
    $cnt++;
  }

  return @jks;
}

sub get_jks {
  # input1: data_ref: [[@conf1],[@conf2],...,[@confN]]
  # input2: binsize
  # if you call this subroutine, you can get jackknif sample
  # as the form ([@Conf1],[@Conf2],...,[Conf{'N/binsize'}])
  my $data_ref = shift;
  my $bin_size = shift;
  my @data = @$data_ref;

  my $num;
#  print "@data\n";
  {
    my @cdata = @{$data_ref->[0]};
    $num = $#cdata; # the number of values in each conf
  }
#  print "$num\n";
  my @total = &new_float_array($num);
  foreach my $cdata_ref ( @$data_ref ) {
    @total = &add_array([@total],$cdata_ref);
#    print "AA $num\n";
  }

  my @jks = ();
  my $flag = 0;
  my $cnt = 0;
  while ( $flag == 0 ) {
    my @tmp = @total;
    my $ctr = 0;
    for ( my $ic = $cnt * $bin_size ; $ic < ($cnt+1) * $bin_size ; $ic++ ) {
      @tmp = &sub_array([@tmp],$data_ref->[$ic]);
      $ctr++;
      if ( $ic == $#data ) {
	$flag = 1;
	last;
      }
    }
    @tmp = &divide_array([@tmp],$#data+1-$ctr);
    push @jks, [@tmp];
    $cnt++;
  }

  return @jks;
}

sub get_jks_multi {
  # input1: data_ref: [[@conf1],[@conf2],...,[@confN]]
  #         @conf{i} = ( [@val{i1}], [@val{i2}],...,[@val{iN}] );
  #         @val{ip} stores some kind of values of interest
  # input2: binsize
  # if you call this subroutine, you can get jackknif sample
  # as the form ([@Conf1],[@Conf2],...,[Conf{'N/binsize'}])
  #   format of [@Conf{i}] is same as [@conf{i}]
  my $data_ref = shift;
  my $bin_size = shift;
  my @data = @$data_ref;
  my $nump;
  my $numv;
  {
    my @cdata = @{$data_ref->[0]};
    $nump = $#cdata; # the number of values in each conf
    my @pdata = @{$cdata[0]};
    $numv = $#pdata;
  }
  my @init = &new_float_array($numv);
  my @total = &new_ref_array($nump,[@init]);

  foreach my $cdata_ref(@data) {
    @total = &add_mtx([@total],$cdata_ref);
  }

  my @jks = ();
  my $flag = 0;
  my $ib = 0;
  while( $flag == 0 ) {
    my @tmp = @total;
    my $ctr = 0;
    for ( my $ic = $ib * $bin_size ; $ic < ($ib+1) * $bin_size ; $ic++ ) {
      @tmp = &sub_mtx([@tmp],$data[$ic]);
      $ctr++;
      if ( $ic == $#data ) {
	$flag = 1;
	last;
      }
    }
    @tmp = &divide_mtx([@tmp],$#data+1-$ctr);
    push @jks, [@tmp];
    $ib++;
  }

  return @jks;
}

sub get_jks_w_err {
  # input1: data_ref: [[@conf1],[@conf2],...,[@confN]]
  # input2: binsize
  # if you call this subroutine, you can get jackknif sample
  # as the form [[@Conf1],[@Conf2],...,[Conf{'N/binsize'}]],[[@Econf1],...,[Econf{'N/binsize'}]]
  my $data_ref = shift;
  my $bin_size = shift;
  my @data = @$data_ref;

  my $num;
  {
    my @cdata = @{$data_ref->[0]};
    $num = $#cdata; # the number of values in each conf
  }
  my @total = &new_float_array($num);
  my @total2 = &new_float_array($num);
  foreach my $cdata_ref ( @$data_ref ) {
    @total = &add_array([@total],$cdata_ref);
    @total2 = &add_array_sq([@total2],$cdata_ref);
  }

  my @jks = ();
  my @jke = ();
  my $flag = 0;
  my $cnt = 0;
  while ( $flag == 0 ) {
    my @tmp = @total;
    my @tmp2 = @total2;
    my $ctr = 0;
    for ( my $ic = $cnt * $bin_size ; $ic < ($cnt+1) * $bin_size ; $ic++ ) {
      @tmp = &sub_array([@tmp],$data_ref->[$ic]);
      @tmp2 = &sub_array_sq([@tmp2],$data_ref->[$ic]);
      $ctr++;
      if ( $ic == $#data ) {
	$flag = 1;
	last;
      }
    }
    @tmp = &divide_array([@tmp],$#data+1-$ctr);
    @tmp2 = &divide_array([@tmp2],$#data-$ctr);

    my @err = &mult_array_each_comp([@tmp],[@tmp]);
    @err = &mult_array([@err],($#data+1-$ctr)/($#data-$ctr));
    @err = &sub_array([@tmp2],[@err]);
    @err = &sqrt_array(@err);

    push @jks, [@tmp];
    push @jke, [@err];
    $cnt++;
  }

  return ( [@jks], [@jke] );
}

sub sort_to_jks_bin {
  # input1 : sort_fname
  # input2 : size of 1conf data in sort data
  # input3 : bin size
  # output : jackknife sample in standard format
  my $sort_fname = shift;
  my $size_cdata = shift;
  my $bs = shift;
  my $nconf = ( -s $sort_fname ) / $size_cdata;
  open (ST,"< $sort_fname") or die "Failed to open $sort_fname";
  binmode ST;
#  print "size of cdata: $size_cdata,  #conf: $nconf\n";
  my @data = ();
  for ( my $ic = 0 ; $ic < $nconf ; $ic++ ) {
    my $buf;
    read (ST, $buf, $size_cdata);
    push @data, [unpack("d*",$buf)];
  }
  close (ST);

  return &Functions::get_jks([@data],$bs);
}

sub save_jks {
  # input1 : file name to save jackknife sample
  # input2 : jackknife sample in the format [[Conf1],[Conf2],...,[ConfNb]]
  # input3 : label for one data [Confi][<label>]
  my $fname = shift;
  my $jks_ref = shift;
  my $label_ref = shift;
  my @jks = @$jks_ref;
  my @label = @$label_ref;
  open (JKSS,"> $fname") or die "Failed to open $fname\n";
  for ( my $ic = 0 ; $ic <= $#jks ; $ic++ ) {
    my @jk = @{$jks[$ic]};
    print "Warning : sub save_jks. $#jk  $#label\n" unless ( $#jk == $#label );
    for ( my $ix = 0 ; $ix <= $#jk ; $ix++ ) {
      print JKSS "$label[$ix] $jk[$ix]\n";
    }
    print JKSS "#---------------------------------------\n";
  }
  print " Jackknife sample saved on $fname in txt format\n";
  close(JKSS);
  return;
}

sub save_jks_multi_nolabel {
  # input1 : file name to save jackknife sample
  # input2 : jackknife sample in the format [Conf1],[Conf2],...,[ConfNb]
  my $fname = shift;
  my @jks = @_;# $jks[$ic]->[$ix]->[$ich]
  open (JKSS,"> $fname") or die "Failed to open $fname\n";
  for ( my $ic = 0 ; $ic <= $#jks ; $ic++ ) {
    my @jk = @{$jks[$ic]};
    for ( my $ix = 0 ; $ix <= $#jk ; $ix++ ) {
      my @jkx = @{$jk[$ix]};
      print JKSS "@jkx\n";
    }
    print JKSS "#---------------------------------------\n";
  }
  print " Jackknife sample saved on $fname in txt format\n";
  close(JKSS);
  return;
}

sub save_jks_bin {
  my $jks_fname = shift;
  my @jks = @_;
  open (JK,"> $jks_fname\n") or die "Failed to open $jks_fname\n";
  binmode JK;
  foreach my $cdata_ref (@jks) {
    print JK pack("d*",@$cdata_ref);
  }
  close (JK);
  print "Jackknife sample saved on $jks_fname\n";
  return;
}

sub read_jks_bin {
  my $jks_fname = shift;
  my $size_cdata = shift;
  my $nconf = ( -s $jks_fname ) / $size_cdata;
  open (JK,"< $jks_fname") or die "Failed to open $jks_fname\n";
  binmode JK;
  my @jks = ();
  for ( my $ic = 0 ; $ic < $nconf ; $ic++ ) {
    my $buf;
    read (JK, $buf, $size_cdata);
    push @jks, [unpack("d*",$buf)];
  }
  close (JK);
  return @jks;
}

sub save_jks_multi_bin {
  my $jks_fname = shift;
  my @jks = @_;
  open (JK,"> $jks_fname\n") or die "Failed to open $jks_fname\n";
  binmode JK;
  foreach my $cdata_ref (@jks) {
    foreach my $pdata_ref (@$cdata_ref) {
      print JK pack("d*",@$pdata_ref);
    }
  }
  close (JK);
  print "Multiple Jackknife sample saved on $jks_fname\n";
  return;
}

sub read_jks_multi_bin {
  my $jks_fname = shift;
  my $Np = shift;
  my $size_pdata = shift;
  my $nconf = ( -s $jks_fname ) / $Np / $size_pdata;
  open (JK,"< $jks_fname") or die "Failed to open $jks_fname\n";
  binmode JK;
  my @jks = ();
  for ( my $ic = 0 ; $ic < $nconf ; $ic++ ) {
    my @cdata = ();
    for ( my $ip = 0 ; $ip < $Np ; $ip++ ) {
      my $buf;
      read (JK, $buf, $size_pdata);
      push @cdata, [unpack("d*",$buf)];
    }
    push @jks, [@cdata];
  }
  close (JK);
  return @jks;
}

sub analyze_jks_single {
  # input  : jackknife sample in the format ($Conf1,$Conf2,...,$ConfNb)
  # output : average & error in the format ($ave,$err)
  my @jks = @_;

  my $ave  = 0;
  my $ave2 = 0;
  foreach my $cdata ( @jks ) {
    $ave  += $cdata;
    $ave2 += $cdata*$cdata;
  }
  $ave  /= $#jks+1;
  $ave2 /= $#jks+1;

  my $err = sqrt($#jks * ($ave2 - $ave*$ave));

  return ( $ave, $err );
}

sub analyze_jks {
  # input  : jackknife sample in the format ([Conf1],[Conf2],...,[ConfNb])
  # output : average & error in the format ([$ave],[$err])
  my @jks = @_;
  my $num;
  {
    my @cdata = @{$jks[0]};
    $num = $#cdata; # number of values in each Conf
  }
  my @ave  = &new_float_array($num);
  my @ave2 = &new_float_array($num);
#  print "AAA\n";
#  my $ctr = 0;
  foreach my $cdata_ref ( @jks ) {
    my @test = @$cdata_ref;
#    $ctr++;
#    print "$ctr $#test\n";
    @ave  = &add_array([@ave],$cdata_ref);
    @ave2 = &add_array_sq([@ave2],$cdata_ref);
  }
#  print "BBB\n";
  @ave  = &divide_array([@ave], $#jks+1);
  @ave2 = &divide_array([@ave2],$#jks+1);

  my @err = &new_float_array($num);
  for ( my $id = 0 ; $id <= $num ; $id++ ) {
    my $test = $#jks * ($ave2[$id] - $ave[$id]*$ave[$id]);
#    if ( abs($test) < $ave[$id] * 1e-15 ) {
    if ( $test < 0 ) {
      $err[$id] = 0;
      next;
    }
#    print "$id $ave[$id] $test\n";
    $err[$id] = sqrt($test);
  }

  return ( [@ave], [@err] );
}

sub analyze_jks_multi {
  # input  : jackknife sample in the format ([Conf1],[Conf2],...,[ConfNb])
  #         @Conf{i} = ( [@val{i1}], [@val{i2}],...,[@val{iN}] );
  #         @val{ip} stores some kind of values of interest
  # output : average & error in the format ([$ave],[$err])
  my @jks = @_;
  my $nump;
  my $numv;
  {
    my @cdata = @{$jks[0]};
    $nump = $#cdata; # number of points in each Confs
    my @pdata = @{$cdata[0]};
    $numv = $#pdata; # number of values in each points
  }
  my @init = &new_float_array($numv);
  my @ave  = &new_ref_array($nump,[@init]);
  my @ave2 = @ave;
  foreach my $cdata_ref ( @jks ) {
    @ave  = &add_mtx([@ave],$cdata_ref);
    @ave2 = &add_mtx_sq([@ave2],$cdata_ref);
  }
  @ave  = &divide_mtx([@ave], $#jks+1);
  @ave2 = &divide_mtx([@ave2],$#jks+1);

  my @err = ();
  for ( my $ip = 0 ; $ip <= $nump ; $ip++ ) {
    my @tmp = ();
    for ( my $iv = 0 ; $iv <= $numv ; $iv++ ) {
      my $sigma = $#jks * ($ave2[$ip]->[$iv] - $ave[$ip]->[$iv]**2);
      unless ($ave[$ip]->[$iv]==0) {
        if ( abs( $sigma/$ave[$ip]->[$iv]**2 ) < 1e-12 ) {
          push(@tmp,0);
          next;
        }
      }
      push(@tmp,sqrt( $sigma ));
    }
    push @err, [@tmp];
  }

  return ( [@ave], [@err], $numv );
}

sub get_cov_mtx_jks {
  # input  : jackknife sample in the format ([Conf1],[Conf2],...,[ConfNb])
  # output : covariance mtx in the format $sigma[p1]->[p2]
  # $sigma[p1]->[p2] = (N_jk - 1)/N_jk
  #          * Σ_ic ( $cdata{ic}[p1] - $ave[p1] ) ( $cdata{ic}[p2] - $ave[p2] )
  my @jks = @_;
  my $num;
  {
    my @cdata = @{$jks[0]};
    $num = $#cdata; # number of values in each Conf
  }
  my @ave  = &new_float_array($num);
  foreach my $cdata (@jks) {
    @ave = &add_array([@ave],$cdata);
  }
  @ave = &divide_array([@ave],$#jks+1);

  my @sigma = @{&new_float_mtx($num,$num)};
  for ( my $ic = 0 ; $ic <= $#jks ; $ic++ ) {
    for ( my $p1 = 0 ; $p1 <= $num ; $p1++ ) {
      for ( my $p2 = 0 ; $p2 <= $num ; $p2++ ) {
        $sigma[$p1]->[$p2]
            += ($jks[$ic]->[$p1] - $ave[$p1]) * ($jks[$ic]->[$p2] - $ave[$p2]);
      }
    }
  }
  @sigma = &mult_mtx([@sigma],$#jks/($#jks+1));

  return @sigma;
}

sub cov_sjks_single {
  # input  : jackknife sample in the format ([@ens1],[@ens2],...,[@ensNe]);
  #     @ens{n} has N_jk real values.  N_jk: independent of n.
  # output : average & error in the format ([@ave],[@cov])
  #          @ave = ($ave1,$ave2,...,$aveNe);
  #     $cov[e1,e2] = (N_jk - 1)/N_jk
  #            * Σ_i^N_jk ( $ens{e1}[i] - $ave{e1} ) ( $ens{e2}[i] - $ave{e2} )
  my @sjks = @_;

  my @ave  = ();
  foreach my $ens_ref (@sjks) {
    my @ens = @$ens_ref;
    my $val = 0;
    foreach my $cval (@ens) {
      $val += $cval;
    }
    $val /= $#ens+1;
    push(@ave,$val);
  }
  # $ave[$ie]

  my @cov = ();
  for ( my $ie = 0 ; $ie <= $#ave ; $ie++ ) {
    my @tmp = ();
    my @cval_i = @{$sjks[$ie]};
    my $ave_i = $ave[$ie];
    for ( my $je = 0 ; $je <= $#ave ; $je++ ) {
      my $val = 0;
      my @cval_j = @{$sjks[$je]};
      my $ave_j = $ave[$je];
      for ( my $i = 0 ; $i <= $#cval_i ; $i++ ) {
        $val += ($cval_i[$i] - $ave_i) * ($cval_j[$i] - $ave_j);
      }
      $val *= $#cval_i / ($#cval_i+1);
      push(@tmp,$val);
    }
    push @cov, [@tmp];
  }
  
  return ( [@ave], [@cov] );
}

sub cov_to_err {
  my @cov = @_;
  my @err = ();
  for ( my $i = 0 ; $i <= $#cov ; $i++ ) {
    push(@err,sqrt($cov[$i]->[$i]));
  }
  print "Error from Cov_Mtx. = @err\n";
  return @err;
}

sub analyze_data {
  # input  : data in the format ([conf1],[conf2],...,[confN]), not jks
  # output : average & error in the format ([$ave],[$err])
  my @data = @_;
  my $num;
  {
    my @cdata = @{$data[0]};
    $num = $#cdata; # number of valuew in each conf
  }
  my @ave  = &new_float_array($num);
  my @ave2 = &new_float_array($num);
  foreach my $cdata_ref ( @data ) {
    @ave  = &add_array([@ave],$cdata_ref);
    @ave2 = &add_array_sq([@ave2],$cdata_ref);
  }
  @ave  = &divide_array([@ave], $#data+1);
  @ave2 = &divide_array([@ave2],$#data+1);

  my @err = &new_float_array($num);
  for ( my $id = 0 ; $id <= $num ; $id++ ) {
    $err[$id] = sqrt( ($ave2[$id] - $ave[$id]*$ave[$id]) / $#data )
	if ( $#data != 0 );
  }

  return ( [@ave], [@err] );
}

sub analyze_data_multi {
  # input  : raw sample in the format ([Conf1],[Conf2],...,[ConfNb])
  #         @Conf{i} = ( [@val{i1}], [@val{i2}],...,[@val{iN}] );
  #         @val{ip} stores some kind of values of interest
  # output : average & error in the format ([$ave],[$err])
  my @jks = @_;
  my $nump;
  my $numv;
  {
    my @cdata = @{$jks[0]};
    $nump = $#cdata; # number of points in each Confs
    my @pdata = @{$cdata[0]};
    $numv = $#pdata; # number of values in each points
  }
  my @init = &new_float_array($numv);
  my @ave  = &new_ref_array($nump,[@init]);
  my @ave2 = @ave;
  foreach my $cdata_ref ( @jks ) {
    @ave  = &add_mtx([@ave],$cdata_ref);
    @ave2 = &add_mtx_sq([@ave2],$cdata_ref);
  }
  @ave  = &divide_mtx([@ave], $#jks+1);
  @ave2 = &divide_mtx([@ave2],$#jks+1);

  my @err = ();
  for ( my $ip = 0 ; $ip <= $nump ; $ip++ ) {
    my @tmp = ();
    for ( my $iv = 0 ; $iv <= $numv ; $iv++ ) {
      my $sigma = ($ave2[$ip]->[$iv] - $ave[$ip]->[$iv]**2);
      unless ($ave[$ip]->[$iv]==0) {
        if ( abs( $sigma/$ave[$ip]->[$iv]**2 ) < 1e-12 ) {
          push(@tmp,0);
          next;
        }
      }
      push(@tmp,sqrt( $sigma ));
    }
    push @err, [@tmp];
  }

  return ( [@ave], [@err], $numv );
}

sub save_result {
  my $result_fname = shift;
  my $ave = shift;
  my $err = shift;
  my @p = @_;

  open (RES,"> $result_fname") or die "Failed to open $result_fname\n";
  for ( my $ip = 0 ; $ip <= $#p ; $ip++ ) {
    print RES sprintf("%s  %e  %e\n",$p[$ip],$ave->[$ip],$err->[$ip]);
  }
  close (RES);
  print "Result written on $result_fname\n";
  return;
}

sub save_result_prec {
  my $result_fname = shift;
  my $ave = shift;
  my $err = shift;
  my @p = @_;

  open (RES,"> $result_fname") or die "Failed to open $result_fname\n";
  for ( my $ip = 0 ; $ip <= $#p ; $ip++ ) {
    print RES sprintf("%s  %.15e  %.15e\n",$p[$ip],$ave->[$ip],$err->[$ip]);
  }
  close (RES);
  print "Result written on $result_fname\n";
  return;
}

sub save_result_multi {
  my $result_fname = shift;
  my $ave = shift;
  my $err = shift;
  my $num_v = shift;
  my @p = @_;

  open (RES,"> $result_fname") or die "Failed to open $result_fname\n";
  for ( my $ip = 0 ; $ip <= $#p ; $ip++ ) {
    print RES sprintf("%s",$p[$ip]);
    for ( my $iv = 0 ; $iv <= $num_v ; $iv++ ) {
      print RES sprintf("  %.8e  %.5e",$ave->[$ip]->[$iv],$err->[$ip]->[$iv]);
    }
    print RES "\n";
  }
  close (RES);
  print "Result written on $result_fname\n";
  return;
}

sub get_rand_bss {
  my $rand_fname = shift;
  my $B = shift;
  my $num_data = shift;

  my @result = ();
  if ( -f $rand_fname ) {
    open(RD,"< $rand_fname") or die "Failed to open $rand_fname\n";
    my @tmp = ();
    while(<RD>) {
      chomp;
      if (/------/) {
	push @result, [@tmp];
	@tmp = ();
	next;
      }
      push(@tmp,$_);
    }
    close(RD);
  } else {
    open(RD,"> $rand_fname") or die "Failed to open $rand_fname\n";
    for ( my $ib = 0 ; $ib < $B ; $ib++ ) {
      my @tmp = ();
      for ( my $ic = 0 ; $ic <= $num_data ; $ic++ ) {
	my $rand = int(rand($num_data+1));
	print RD "$rand\n";
	push(@tmp,$rand);
      }
      print RD "-------------------------\n";
      push @result, [@tmp];
    }
    close(RD);
  }
  return @result;
}

sub get_bss {
  # input1: data_ref: [[@conf1],[@conf2],...,[@confN]]
  # input2: bootstrap size
  # if you call this subroutine, you can get bootstrap sample
  # as the form ([@Conf1],[@Conf2],...,[Conf{'N/binsize'}])
  my $data_ref = shift;
  my $B = shift;
  my $rand_fname = (shift||"bootstrap.rand");
  my @data = @$data_ref;

  my $num;
  {
    my @cdata = @{$data_ref->[0]};
    $num = $#cdata; # the number of values in each conf
  }
  my @init = &new_float_array($num);

  my @rand_nums = &get_rand_bss($rand_fname,$B,$#data);
  my @bss = ();
  for ( my $ib = 0 ; $ib < $B ; $ib++ ) {
    my @bdata = @init;
    my @brands = @{$rand_nums[$ib]};
    for ( my $ic = 0 ; $ic <= $#brands ; $ic++ ) {
      my $nc = $brands[$ic];
      @bdata = &add_array([@bdata],$data[$nc]);
    }
    @bdata = &divide_array([@bdata],$#brands+1);
    push @bss, [@bdata];
  }

  return @bss;
}

sub get_bss_multi {
  # input1: data_ref: [[@conf1],[@conf2],...,[@confN]]
  #         @conf{i} = ( [@val{i1}], [@val{i2}],...,[@val{iN}] );
  #         @val{ip} stores some kind of values of interest
  # input2: Number of bin
  # if you call this subroutine, you can get bootstrap sample
  # as the form ([@Conf1],[@Conf2],...,[Conf{B}])
  #   format of [@Conf{i}] is same as [@conf{i}]
  my $data_ref = shift;
  my $B = shift;
  my $rand_fname = shift;
  my @data = @$data_ref;
  my $nump;
  my $numv;
  {
    my @cdata = @{$data_ref->[0]};
    $nump = $#cdata; # the number of values in each conf
    my @pdata = @{$cdata[0]};
    $numv = $#pdata;
  }
  my @inita = &new_float_array($numv);
  my @init = &new_ref_array($nump,[@inita]);

  my @rand_nums = &get_rand_bss($rand_fname,$B,$#data);
  my @bss = ();
  for ( my $ib = 0 ; $ib < $B ; $ib++ ) {
    my @bdata = @init;
    my @brands = @{$rand_nums[$ib]};
    for ( my $ic = 0 ; $ic <= $#brands ; $ic++ ) {
      my $nc = $brands[$ic];
      @bdata = &add_mtx([@bdata],$data[$nc]);
    }
    @bdata = &divide_mtx([@bdata],$#brands+1);
    push @bss, [@bdata];
  }

  return @bss;
}

sub analyze_bss {
  # input  : bootstrap sample in the format ([Conf1],[Conf2],...,[ConfNb])
  # output : average & error in the format ([$ave],[$err])
  my @bss = @_;
  my $num;
  {
    my @cdata = @{$bss[0]};
    $num = $#cdata; # number of values in each Conf
  }
  my @ave  = &new_float_array($num);
  my @ave2 = &new_float_array($num);
  foreach my $cdata_ref ( @bss ) {
    my @test = @$cdata_ref;
    @ave  = &add_array([@ave],$cdata_ref);
    @ave2 = &add_array_sq([@ave2],$cdata_ref);
  }
  @ave  = &divide_array([@ave], $#bss+1);
  @ave2 = &divide_array([@ave2],$#bss+1);

  my @err = &new_float_array($num);
  for ( my $id = 0 ; $id <= $num ; $id++ ) {
#    print "$id\n";
    if ( abs(($ave2[$id] - $ave[$id]*$ave[$id])/($ave2[$id] + $ave[$id]*$ave[$id])) < 1e-10 ) {
      $err[$id] = 0;
      next;
    }
    $err[$id] = sqrt(($ave2[$id] - $ave[$id]*$ave[$id]));
  }

  return ( [@ave], [@err] );
}

sub analyze_bss_multi {
  # input  : bootstrap sample in the format ([Conf1],[Conf2],...,[ConfNb])
  #         @Conf{i} = ( [@val{i1}], [@val{i2}],...,[@val{iN}] );
  #         @val{ip} stores some kind of values of interest
  # output : average & error in the format ([$ave],[$err])
  my @bss = @_;
  my $nump;
  my $numv;
  {
    my @cdata = @{$bss[0]};
    $nump = $#cdata; # number of points in each Confs
    my @pdata = @{$cdata[0]};
    $numv = $#pdata; # number of values in each points
  }
  my @init = &new_float_array($numv);
  my @ave  = &new_ref_array($nump,[@init]);
  my @ave2 = @ave;
  foreach my $cdata_ref ( @bss ) {
    @ave  = &add_mtx([@ave],$cdata_ref);
    @ave2 = &add_mtx_sq([@ave2],$cdata_ref);
  }
  @ave  = &divide_mtx([@ave], $#bss+1);
  @ave2 = &divide_mtx([@ave2],$#bss+1);

  my @err = ();
  for ( my $ip = 0 ; $ip <= $nump ; $ip++ ) {
    my @tmp = ();
    for ( my $iv = 0 ; $iv <= $numv ; $iv++ ) {
      my $sqrt = $ave2[$ip]->[$iv] - $ave[$ip]->[$iv]**2;
      if ( abs($sqrt) <= 1e-12*$ave2[$ip]->[$iv] && $sqrt < 0 ) {
        $sqrt = 0;
      } elsif ( $sqrt < 0 ) {
#        print ($ave2[$ip]->[$iv] - $ave[$ip]->[$iv]**2) . " $sqrt\n";
        die;
      } else {
        $sqrt = sqrt($sqrt);
      }
      push(@tmp,$sqrt);
    }
    push @err, [@tmp];
  }

  return ( [@ave], [@err], $numv );
}

sub round_off {
  my $val = shift;
  my $digit = shift;
  # in this case, $digit is last digit which you want to keep
  # ex) if you want to map 0.1893744 to 0.189, you should input $digit as -3
  #                       -34.578315 to -34.6, you should input $digit as -1
  my $sign = $val / abs($val);

  my $a = $val*10**(-$digit)+0.5*$sign;

  return (int($val*10**(-$digit)+0.5*$sign))* 10 **($digit);
}

sub round_off_ndigit {
  my $val = shift;
  my $digit = shift;
  # in this case, $digit is # of digits which you want to keep
  # ex) if you want to map 0.1893744 to 0.189, you should input $digit = 3
  #                       -34.578315 to -34.6, you should input $digit = 3
  #                       -34.578315 to -35, you should input $digit = 2

  my $sign  = $val / abs($val);
  my $min_odr = &gauss_bracket(&log10(abs($val))) - $digit + 1;
  return &round_off($val,$min_odr);
}

sub tablize {
  my $val = shift;
  my $err = shift;
  my $err_digit = ( shift || 2 );

#  print "tablize $val & $err\n";
  my $odr_val = &gauss_bracket(&log10(abs($val)));
  my $odr_err = &gauss_bracket(&log10($err));
  return ("-","-") if ( $odr_val < $odr_err );

  my $digit_decimal = $odr_err - $err_digit + 1;
  my $all_digit = 2 - $digit_decimal;# 2 is because of "0" & "."
  $all_digit += $odr_val if ( $odr_val > 0 );

  my $valr = sprintf("%*.*f", $all_digit, - $digit_decimal,
		                   &round_off($val,$digit_decimal));
  $valr = sprintf("-%*.*f", $all_digit, - $digit_decimal,
		               abs(&round_off($val,$digit_decimal)))
                             if ( $val < 0 );
  my $errr = &round_off($err,$digit_decimal) * 10 **(-$digit_decimal);
  my @result = ($valr,$errr);

  return @result;
}

sub get_texf_number {
  my $val = shift;
  my $err = shift;
  my $err_digit = (shift||2);
  # Format: $*.**(**)\times10^{*}$

  my $digit = 1;
  $digit = &gauss_bracket(&log10(abs($val))) unless ($val==0);
  my $sval = $val / 10**$digit;

  my $edigit = &gauss_bracket(&log10(abs($err)));
  if ( $val == 0 ) {
    my $result = sprintf("\$0(%.f)\\times10^\{%d\}\$",($err/10**$edigit),$edigit);
    return $result;
  }
  my $eff_digit = $digit - $edigit + $err_digit;
  if ( $eff_digit < 1 ) {
    my $tval = sprintf("%.f",$sval);
    my $terr = sprintf("%.f",($err/10**$digit));
    my $result = sprintf("\$%s(%s)\\times10^\{%d\}\$",$tval,$terr,$digit);
    $result = sprintf("\$%s(%s)\$",$tval,$terr) if ( $digit == 0 );
    return $result;
  }
  my $tmp = "1." . ($eff_digit-1);
  my $tval = sprintf("%${tmp}f",$sval);
  my $terr = sprintf("%.f",$err*10**($eff_digit-$digit-1));
  if ( $digit == $edigit ) {
    $terr = sprintf("%${tmp}f",($err/10**$edigit));
  }
  my $result = sprintf("\$%s(%s)\\times10^\{%d\}\$",$tval,$terr,$digit);
  $result = sprintf("\$%s(%s)\$",$tval,$terr) if ( $digit == 0 );
  return $result;
}

sub get_gpf_number_woerr {
  my $val = shift;
  my $eff_digit = (shift||3);
  # Format: *.** {/Symbol \264} 10^{*}

  my $digit = 1;
  $digit = &gauss_bracket(&log10(abs($val))) unless ($val==0);
  my $sval = $val / 10**$digit;

  if ( $val == 0 ) {
    my $tmp = "0" . ($eff_digit-1);
    my $result = sprintf("0.%${tmp}d",0);
    return $result;
  }
  my $tmp = "1." . ($eff_digit-1);
  my $tval = sprintf("%${tmp}f",$sval);
  my $result = sprintf("%s \{/Symbol \\264\} 10^\{%d\}",$tval,$digit);
  $result = sprintf("%s",$tval) if ( $digit == 0 );
  return $result;
}

##
# Taylor Expansion
##
our $Func_exp = "exp";
sub get_func_exp { return $Func_exp };
sub taylor_expand_exp {
  # get coefficients of taylor expansion of the function
  # exp(x)
  my $deg = ( shift || 10 );
  my @result = (1);
  my $tmp = 1;
  for ( my $i = 1 ; $i <= $deg ; $i++ ) {
    $tmp /= $i;
    push(@result,$tmp);
  }
  return @result;
}

our $Func_1_1px = "1/(1+x)";
sub get_func_1_1px { return $Func_1_1px };
sub taylor_expand_1_1px {
  # get coefficients of taylor expansion of the function
  # 1 / ( 1 + x )
  my $deg = ( shift || 10 );
  my @result = ();

  for ( my $i = 0 ; $i <= $deg ; $i++ ) {
    push(@result,(-1)**$i);
  }
  return @result;
}

our $Func_pow_1px_a = "(1+x)^a";
sub get_func_pow_1px_a { return $Func_pow_1px_a };
sub taylor_expand_pow_1px_a {
  # get coefficients of taylor expansion of the function
  # ( 1 + x )^a
  my $deg = shift;
  my $a = shift;
  my @result = ();
  if ( $a == int($a) && $a >= 0 ) {
    for ( my $i = 0 ; $i <= &min($a,$deg) ; $i++ ) {
      push(@result,&factorial($a)/&factorial($i)/&factorial($a-$i));
    }
    while( $#result < $deg ) {
      push(@result,0.0);
    }
    return @result;
  }

  my $tmp = 1;
  push (@result,1);
  for ( my $i = 1 ; $i <= $deg ; $i++ ) {
    $tmp *= ($a+1-$i);
    $tmp /= $i;
    push(@result,$tmp);
  }
  return @result;
}

our $Func_ln1px = "ln(1+x)";
sub get_func_ln1px { return $Func_ln1px };
sub taylor_expand_ln1px {
  # get coefficients of taylor expansion of the function
  # ln ( 1 + x )
  my $deg = ( shift || 10 );
  my @result = (0);
  for ( my $i = 1 ; $i <= $deg ; $i++ ) {
    push(@result,1/$i*(-1)**($i+1));
  }
  return @result;
}

our $Func_pow_x_n = "x^n";
sub get_func_pow_x_n { return $Func_pow_x_n }
sub taylor_expand_pow_x_n {
  my $deg = shift;
  my $n = shift;
  my @result = &new_int_array($deg);
  $result[$n] = 1;

  return @result;
}

sub taylor_expansion {
  my $function = shift;
  my $deg = shift;
  my $c1 = ( shift || 0 );
  return &taylor_expand_exp($deg) if ($function eq $Func_exp);
  return &taylor_expand_1_1px($deg) if ($function eq $Func_1_1px);
  return &taylor_expand_pow_x_n($deg,$c1) if ($function eq $Func_pow_x_n);
  return &taylor_expand_pow_1px_a($deg,$c1) if ($function eq $Func_pow_1px_a);
  return &taylor_expand_ln1px($deg) if ($function eq $Func_ln1px);
  die;
}

our $Func_polynomial = "polynomial";
sub get_func_polynomial { return $Func_polynomial }

sub get_add_combination {
  my $target_value = shift;
  my $num_divide = shift;

  my @result = ();
  my $pos_1dim = $target_value - $num_divide + 1;
  my $total_pattern = $pos_1dim ** $num_divide;
  for ( my $i = 0 ; $i < $total_pattern ; $i++ ) {
    my $tmp1;
    my $tmp2 = $i;
    my @cand = ();
    my $sum = 0;
    for ( my $idiv = 1 ; $idiv <= $num_divide ; $idiv++ ) {
      $tmp1 = $tmp2 % $pos_1dim + 1;
      $tmp2 = int( $tmp2 / $pos_1dim );
      $sum += $tmp1;
      push (@cand,$tmp1)
    }
    next unless ( $sum == $target_value );
#p    print "@cand\n";
    push @result, [@cand];
  }
  return @result;
}

sub order_taylor_expansion {
  # This subroutine calculates a taylor expansion of any polinomials
  # and 
  # get 
  # coefficients of consequence polinomials
  # e.g.
  # P(x) = a0 + a1 x + a2 x^2 + ...
  # P'(x) = 1 + b1 x + b2 x^2 + ..., bi = ai/a0
  # P~(x) = P'(x) - 1
  # 1/P(x) = (1/a0) ( C0 + C1 P~(x) + C2 P~(x)^2 + ... )
  #        = R0 + R1 x + R2 x^2 + ...
  # This subroutine returns ( R0, R1, R2, ... )
  my $in_poly_ref = shift;
  my $type = shift;
  my $deg;
  {
    my @a = @$in_poly_ref;
    $deg = ( shift || $#a );
  }
  my $c1 = ( shift || 0 );

  my $poly_ref = $in_poly_ref;
  return @$poly_ref if ( $type eq $Func_polynomial );
  my $overall_factor = 0.0;
  if ( $type eq $Func_exp ) {
    $overall_factor = exp($poly_ref->[0]);
  } elsif ( $type eq $Func_1_1px ) {
    $overall_factor = 1/$poly_ref->[0];
    my @poly = &divide_array($poly_ref,$poly_ref->[0]);
    $poly_ref = [@poly];
  } elsif ( $type eq $Func_ln1px ) {
    $overall_factor = log($poly_ref->[0]);
    my @poly = &divide_array($poly_ref,$poly_ref->[0]);
    $poly_ref = [@poly];
  } elsif ( $type eq $Func_pow_1px_a ) {
    if ( $poly_ref->[0] == 0 ) {
	print "You are going to perform an expansion at singular point\n"
	  unless ( $c1 == int($c1) );
      $type = $Func_pow_x_n;
    } else {
      $overall_factor = $poly_ref->[0] **$c1;
      my @poly = &divide_array($poly_ref,$poly_ref->[0]);
      $poly_ref = [@poly];
    }
  } else {
    print "You must add in sub ``order_taylor_expansion'' about function type: "
	. "$type\n" unless ( $type eq $Func_pow_x_n );
  }

  my @taylor = &taylor_expansion($type,$deg,$c1);
##  print "$#taylor\n";
#p  print "BBB @taylor\n" if ( $type eq $Func_pow_x_n );
  my @result = ($taylor[0]);
  for ( my $ideg = 1 ; $ideg <= $deg ; $ideg++ ) {
    my $tmp = 0.0;
    for ( my $i = 1 ; $i <= $ideg ; $i++ ) {
      my @comb = &get_add_combination($ideg,$i);
      foreach my $combi (@comb) {
	my $tmp2 = 1.0;
#p	print "ideg: $ideg, it: $i, " if ( $type eq $Func_pow_x_n );
	foreach my $idx (@$combi) {
#p	  print $idx . " ";
	  $tmp2 *= $poly_ref->[$idx];
	}
#p	print "\n";
	$tmp += $tmp2 * $taylor[$i];
      }
    }
    push (@result,$tmp);
  }
  if ( $type eq $Func_ln1px ) {
    for ( my $i = 0 ; $i <= $#result ; $i++ ) {
      $result[$i] += $overall_factor;
    }
  } elsif ( $type eq $Func_pow_x_n ) {
  } else {
    @result = &mult_array([@result],$overall_factor);
  }
  return @result;
}

sub polynomial_product {
  # input1: ref of array of coefficients of polynomial1
  # input2: ref of array of coefficients of polynomial2
  my $poly1_ref = shift;
  my $poly2_ref = shift;
  my $deg;
  {
    my @a1 = @$poly1_ref;
    my @a2 = @$poly2_ref;
#p    print " a1: @a1\n a2: @a2\n";
    print "Warning: sub polynomial_product\n" unless ( $#a1 == $#a2 );
    $deg = ( shift || &min($#a1,$#a2) );
  }

  my @result = ();
  for ( my $ideg = 0 ; $ideg <= $deg ; $ideg++ ) {
    my @comb = &get_add_combination($ideg+2,2);
    my $tmp = 0.0;
    foreach my $combi (@comb) {
      $tmp += $poly1_ref->[$combi->[0]-1] * $poly2_ref->[$combi->[1]-1];
#      $tmp += $poly1_ref->[$test[0]-1] * $poly2_ref->[$test[1]-1];
    }
    push (@result,$tmp);
  }
  return @result;
}

sub polynomial_calculation {
  # This subroutine calculates the product of 2 functions which are
  # polynomial or able to be performed taylor expansion
  #
  # input1: function1   input4: type of function1
  # input2: function2   input5: type of function2
  # input3: max degree
  my $func1_ref = shift;
  my $func2_ref = shift;
  my $deg;
  {
    my @a1 = @$func1_ref;
    my @a2 = @$func2_ref;
    print "Warning: sub polynomial_calculation\n" unless ( $#a1 == $#a2 );
    $deg = ( shift || &min($#a1,$#a2) );
  }
  my $type1 = ( shift || $Func_polynomial );
  my $type2 = ( shift || $Func_polynomial );

  my @func1 = &order_taylor_expansion($func1_ref,$type1,$deg);
  my @func2 = &order_taylor_expansion($func2_ref,$type2,$deg);

  return &polynomial_product([@func1],[@func2],$deg);
}

sub polynomial_integration {
  my $poly_ref = shift;
  my $int_const = ( shift || 0 );
  my @poly = @$poly_ref;
  my $deg = ( shift || $#poly );

  my @result = ($int_const);
  for ( my $i = 0 ; $i <= $deg ; $i++ ) {
    push (@result,$poly[$i]/($i+1));
  }
  return @result;
}

sub compose_poly {
  # If you call this subroutine with two input polynomials P1 & P2,
  # you can get the new polynomials as a composite function P1(P2(x))
  my $P1_ref = shift;
  my $P2_ref = shift;
  my $deg = shift;

  my @P1 = @$P1_ref;
#  print "P1: @P1\n";
  my @P2 = @$P2_ref;
#  print "P2: @P2\n";
  my @result = &new_float_array($deg);
  $result[0] = $P1_ref->[0];
#  print "AAA @result\n";
  for ( my $i = 1 ; $i <= $deg ; $i++ ) {
    my @Q = &order_taylor_expansion($P2_ref,$Func_pow_x_n,$deg,$i);
#    print "BBB $i @Q\n";
    @Q = &mult_array([@Q],$P1_ref->[$i]);
#    print "CCC $i @Q\n";
    @result = &add_array([@result],[@Q]);
#    print "DDD $i @Q\n";
  }
#  print "ZZZ @result\n";
  return @result;
}

sub compose_poly_test {
  # If you call this subroutine with two input polynomials P1 & P2,
  # you can get the new polynomials as a composite function P1(P2(x))
  my $P1_ref = shift;
  my $P2_ref = shift;
  my $deg = shift;

  my @P1 = @$P1_ref;
  print "P1: @P1\n";
  my @P2 = @$P2_ref;
  print "P2: @P2\n";
  my @result = &new_float_array($deg);
  $result[0] = $P1_ref->[0];
  print "AAA @result\n";
  for ( my $i = 1 ; $i <= $deg ; $i++ ) {
    my @Q = &order_taylor_expansion($P2_ref,$Func_pow_x_n,$deg,$i);
    print "BBB $i @Q\n";
    @Q = &mult_array([@Q],$P1_ref->[$i]);
    print "CCC $i @Q\n";
    @result = &add_array([@result],[@Q]);
    print "DDD $i @Q\n";
  }
  print "ZZZ @result\n";
  return @result;
}

#sub calc_prime {
#  my $func = shift;
#  while( $func =~ /int()/ ) {
#
#  }
#}

sub calculate {
  my $func = shift;
  # This subroutine calculate +,-,*,/,^(**),e
  # e.g. input 
  my $tmp = $func;
  my @vals = ();
  my @oprs = ();
  while( $tmp =~ /^([\+\-]*[^\+\-\*\/\^e]+)([\+\-\*\/\^e])(.+)$/ ) {
    push(@vals,$1);
    push(@oprs,$2);
    $tmp = $3;
  }
  push(@vals,$tmp);

  my $io = &inverse_array_char("e",[@oprs]);
  while( $io >= 0 ) {
    my @valst = ();
    my @oprst = ();
    for ( my $i = 0 ; $i <= $#oprs ; $i++ ) {
      push(@valst,$vals[$i]) if ( $i < $io );
      push(@valst,$vals[$i]."e".$vals[$i+1]) if ( $i == $io );
      push(@valst,$vals[$i+1]) if ( $i > $io );
      push(@oprst,$oprs[$i]) unless ( $i == $io );
    }
    @vals = @valst;
    @oprs = @oprst;
    $io = &inverse_array_char("e",[@oprs]);
  }
  $io = &inverse_array_char("^",[@oprs]);
  while( $io >= 0 ) {
    my @valst = ();
    my @oprst = ();
    for ( my $i = 0 ; $i <= $#oprs ; $i++ ) {
      push(@valst,$vals[$i]) if ( $i < $io );
      push(@valst,$vals[$i]**$vals[$i+1]) if ( $i == $io );
      push(@valst,$vals[$i+1]) if ( $i > $io );
      push(@oprst,$oprs[$i]) unless ( $i == $io );
    }
    @vals = @valst;
    @oprs = @oprst;
    $io = &inverse_array_char("^",[@oprs]);
  }
  $io = &inverse_array_char("*",[@oprs]);
  while( $io >= 0 ) {
    my @valst = ();
    my @oprst = ();
    for ( my $i = 0 ; $i <= $#oprs ; $i++ ) {
      push(@valst,$vals[$i]) if ( $i < $io );
      push(@valst,$vals[$i]*$vals[$i+1]) if ( $i == $io );
      push(@valst,$vals[$i+1]) if ( $i > $io );
      push(@oprst,$oprs[$i]) unless ( $i == $io );
    }
    @vals = @valst;
    @oprs = @oprst;
    $io = &inverse_array_char("*",[@oprs]);
  }
  $io = &inverse_array_char("/",[@oprs]);
  while( $io >= 0 ) {
    my @valst = ();
    my @oprst = ();
    for ( my $i = 0 ; $i <= $#oprs ; $i++ ) {
      push(@valst,$vals[$i]) if ( $i < $io );
      push(@valst,$vals[$i]/$vals[$i+1]) if ( $i == $io );
      push(@valst,$vals[$i+1]) if ( $i > $io );
      push(@oprst,$oprs[$i]) unless ( $i == $io );
    }
    @vals = @valst;
    @oprs = @oprst;
    $io = &inverse_array_char("/",[@oprs]);
  }
  $io = &inverse_array_char("-",[@oprs]);
  while( $io >= 0 ) {
    my @valst = ();
    my @oprst = ();
    for ( my $i = 0 ; $i <= $#oprs ; $i++ ) {
      push(@valst,$vals[$i]) if ( $i < $io );
      push(@valst,$vals[$i]-$vals[$i+1]) if ( $i == $io );
      push(@valst,$vals[$i+1]) if ( $i > $io );
      push(@oprst,$oprs[$i]) unless ( $i == $io );
    }
    @vals = @valst;
    @oprs = @oprst;
    $io = &inverse_array_char("-",[@oprs]);
  }
  $io = &inverse_array_char("+",[@oprs]);
  while( $io >= 0 ) {
    my @valst = ();
    my @oprst = ();
    for ( my $i = 0 ; $i <= $#oprs ; $i++ ) {
      push(@valst,$vals[$i]) if ( $i < $io );
      push(@valst,$vals[$i]+$vals[$i+1]) if ( $i == $io );
      push(@valst,$vals[$i+1]) if ( $i > $io );
      push(@oprst,$oprs[$i]) unless ( $i == $io );
    }
    @vals = @valst;
    @oprs = @oprst;
    $io = &inverse_array_char("+",[@oprs]);
  }
  return $vals[0];
}

sub enhance_vector {
  # $in[i]->[j]  =>  $out[i*Nj+j]
  my @in = @_;
  my @out = ();
  foreach my $in1 (@in) {
    my @part = @$in1;
    foreach my $val (@part) {
      push(@out,$val);
    }
  }
  return @out;
}

1;
