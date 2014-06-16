#!/usr/bin/perl

use strict;
use warnings;

my $defaultOpt = {
	maxBasesPerRow => 20000,		# Means 1 mitochondrial genome (~16000 bases) will fit on a single row
	xTicks => 4000
};

sub ceil($) {
	my ($x) = @_;

	return $x if $x == int $x;
	return $x + 1;
}

my $psHeader = <<THE_END;
%!PS

/inch { 72 mul } bind def

% A4 settings
/height 11.7 inch def
/width 8.3 inch def

% Decent margins
/hmargin 1 inch def
/vmargin 1 inch def

/top height vmargin sub def
/left hmargin def

/rowheight 2 inch def
/rowwidth width hmargin 2 mul sub def
/hplotmargin 0.25 inch def		% Room for tick marks and y-axis labels, also used on right end to allow room for rightmost x-axis label
/vplotmargin 0.25 inch def		% Room for tick marks and x-axis labels at bottom
/titleheight 0.5 inch def
/maxbasesperrow 20000 def
/ticklen 0.1 inch def
/axislabelfontsize 10 def
/axislabelfont {
%	/Times-Roman findfont axislabelfontsize scalefont setfont
	/Helvetica findfont axislabelfontsize scalefont setfont
} bind def

/Times-Roman findfont 12 scalefont setfont

/topgrey 0.75 def
/bottomgrey 0.25 def

/plotwidth
	rowwidth
	hplotmargin 2 mul sub
def

/plotheight
	rowheight
	titleheight sub
	vplotmargin sub
def

THE_END

my $solidBlackRect = <<THE_END;

% The workhorse function: Plot one bar, and take one step to the right.
% Takes 1 parameter, the coverage value.
% NOTE: We can't use bind here, because w and h need to be looked up
% *after* this function is defined.  This will make things slightly
% slower, but who cares.
/p {
	currentpoint
	newpath
	moveto

	h mul
	dup
	0 exch rlineto
	w 0 rlineto
	neg 0 exch rlineto
	closepath
	w 0 rmoveto
	currentpoint
	fill
	moveto
} def
THE_END

my $psFunctions = <<THE_END;
% Move the origin to the origin of the axes in the first row on the page.
/newpage {
	left hplotmargin add
	top titleheight sub plotheight sub
	translate
	0 0 moveto
} bind def

% Move the origin down one row.
/downrow {
	0
	0 rowheight sub
	translate
	0 0 moveto
} bind def

% Takes 1 parameter, the string.
/title {
	0 setgray
%	/Times-Roman findfont 16 scalefont setfont
	/Helvetica findfont 14 scalefont setfont
	0 plotheight ticklen add moveto
	show
	0 0 moveto
} bind def

% Takes 1 parameter, the width in bases, and draws a pair of axes.  The plot always has height plotheight.
/axes {
	0 setgray
	w mul
	0 moveto
	0 0 lineto
	0 plotheight lineto
	stroke
	0 0 moveto
} def

% Draw a tick mark and coverage label on the y axis.
% Takes 1 parameter, the string label.  The height is taken from the current point.  The label will generally "mean" the same thing as the height, e.g. "(42)" and "42",
% but we pass it separately because PostScript doesn't have facilities for dealing with decimal places nicely.
/ytick {
	gsave
	0 setgray
	ticklen neg 0 rlineto
%	ticklen neg axislabelfontsize -0.5 mul rmoveto
	ticklen neg axislabelfontsize -0.3 mul rmoveto		%HACK: "Should" be -0.5, but this makes the text too low because of space at the top...
	currentpoint
	stroke
	moveto
	axislabelfont
%	gsave		%DEBUG
	showrightjustified
%	grestore		%DEBUG
%	0 axislabelfontsize rlineto		%DEBUG
%	currentpoint
%	stroke
%	moveto
%	(XXX123XXX) showrightjustified		%DEBUG
	grestore
} def

% Draws a tick on the y axis.
% Takes 2 parameters, the string label and the coverage value.
/drawytick {
	gsave
	h mul
	0 exch moveto
	ytick
	grestore
} def

% Draws a dashed line at the given height, and adds a tick mark on the y axis.
% Takes 3 parameters:
% - The y axis label (specified explicitly to enable control of decimal places)
% - The width in bases
% - The average coverage value
/drawaverage {
	h mul
	dup

	0 setgray

	gsave
	0 exch moveto
	[3 3] 0 setdash
	exch
	w mul
	0 rlineto
	stroke
	grestore

	0 exch moveto
	ytick
	0 0 moveto
} def

% Display a text string centred at the current point.
% Takes 1 parameter, the string to display.
/showcentred {
	dup
	stringwidth
	pop
	-0.5 mul 0 rmoveto
	show
} bind def

/showrightjustified {
	dup
	stringwidth
	pop
	neg 0 rmoveto
	show
} bind def

% Draw a tick mark and location label on the x axis.
% Takes 1 parameter, the reference genome location to use as a label.  The placement is done using the PostScript current point location.
/xtick {
	gsave
	0 setgray
	w 0.5 mul 0 rmoveto
	0 ticklen neg rlineto
	currentpoint
	stroke
	moveto
	0 ticklen axislabelfontsize add neg rmoveto
	axislabelfont
	showcentred
	grestore
} def

% Convert a raw coverage value to a greyscale value in the range bottomgrey to topgrey.
% Can't use bind because we don't know maxcoverage yet.
/converttogrey {
	topgrey bottomgrey sub
	mul
	maxcoverage div
	bottomgrey add
	1 exch sub
} def

% The workhorse function: Plot one bar, and take one step to the right.
% Takes 1 parameter, the coverage value.
% NOTE: We can't use bind here, because w and h need to be looked up
% *after* this function is defined.  This will make things slightly
% slower, but who cares.
/p {
	currentpoint
	newpath
	moveto

	dup
%	invmaxcoverage mul
%	maxcoverage div
	converttogrey
	setgray

	h mul
	dup
	0 exch rlineto
	w 0 rlineto
	neg 0 exch rlineto
	closepath
	w 0 rmoveto
	currentpoint
	fill
	moveto
} def
THE_END

sub ensureHeaderPrinted($) {
	my ($g) = @_;
	my $outF = $g->{outF};		#HACK: Only needed because Perl's crappy syntax forbids using $g->{outF} as the "indirect" argument to print...

	if (!defined $g->{headerPrinted}) {
		print $outF $psHeader, $psFunctions;
		$g->{headerPrinted} = 1;
	}
}

sub newPage($) {
	my ($g) = @_;
	my $outF = $g->{outF};		#HACK: Only needed because Perl's crappy syntax forbids using $g->{outF} as the "indirect" argument to print...

	ensureHeaderPrinted $g;
	if (!defined $g->{iPage}) {
		# Output the PostScript header and function definitions
		$g->{iPage} = 0;
	} else {
		++$g->{iPage};
		print $outF "showpage\n";
	}

	print $outF "newpage\n";
	$g->{iRow} = 0;
}

#sub newRow($$$) {
#	my ($iBase, $nBases, $g) = @_;
sub newRow($) {
	my ($g) = @_;
	my $outF = $g->{outF};		#HACK: Only needed because Perl's crappy syntax forbids using $g->{outF} as the "indirect" argument to print...

	if (defined $g->{iRow} && $g->{iRow} < $g->{maxRows}) {
		print $outF "downrow\n";
		++$g->{iRow};
	} else {
		newPage $g;
	}

#	print $outF "$nBases axes\n";
}

sub finalise($) {
	my ($g) = @_;

	my $outF = $g->{outF};		#HACK: Only needed because Perl's crappy syntax forbids using $g->{outF} as the "indirect" argument to print...
	print $outF "showpage\n";
}

sub drawRest($$$$$$$) {
	my ($iFirstBase, $nBases, $basesPerRow, $maxCov, $avgCov, $opt, $g) = @_;
	my $outF = $g->{outF};

#	print $outF "$nBases $avgCov drawaverage\n";
	my $prettyAvgCov = sprintf "%.2f", $avgCov;
	print $outF "($prettyAvgCov) $nBases $avgCov drawaverage\n";
	print $outF "($maxCov) $maxCov drawytick\n";		#HACK
	print $outF "$nBases axes\n";

#	for (my $i = int(($iFirstBase - $opt->{startBase}) % $basesPerRow / $opt->{xTicks}) * $opt->{xTicks} + 1; $i < $iFirstBase + $nBases; $i += $opt->{xTicks}) {
#		my $x = $i - 
#		print $outF "
#	}
}

# Safely encode a string for PostScript.
sub psString($) {
	local ($_) = @_;
	s/([()\\])/\\$1/g;
	return $_;
}

# What it says on the tin.  Returns undef if no arguments are passed.
sub max(@) {
	my $max = shift;
	foreach (@_) {
		$max = $_ if $_ > $max;
	}

	return $max;
}

# What it says on the tin.  Returns undef if no arguments are passed.
sub min(@) {
	my $min = shift;
	foreach (@_) {
		$min = $_ if $_ < $min;
	}

	return $min;
}

# $g is the "graphics context", which is updated across calls.
sub plot($$$) {
	my ($inFn, $opt, $g) = @_;
	my $outF = $g->{outF};		#HACK: Only needed because Perl's crappy syntax forbids using $opt->{outF} as the "indirect" argument to print...

	print STDERR "Processing $inFn...\n";

	# Pass 1 determines maximum and average coverage levels
	my $maxCov = 0;
	my $avgCov = 0;
	my $nBasesCovered = 0;
	{
		local $_;
		open my $inF, "<", $inFn or die "Could not open '$inFn' for input: $!";
		my ($loc, $cov);
		while (<$inF>) {
			(undef, $loc, undef, $cov) = split /\t/;
			last if defined $opt->{endBase} && $loc > $opt->{endBase};
			$opt->{startBase} = $loc if !defined $opt->{startBase};
			if ($loc >= $opt->{startBase}) {
				$maxCov = $cov if $cov > $maxCov;
				$avgCov += $cov;
				++$nBasesCovered;
			}
		}
		close $inF;
		$opt->{endBase} = $loc if !defined $opt->{endBase};
	}

	my $len = $opt->{endBase} - $opt->{startBase} + 1;
	my $nRows = ceil($len / $opt->{maxBasesPerRow});

	$avgCov /= $len;
	print STDERR "Reference sequence range: $opt->{startBase}-$opt->{endBase}\n";
	print STDERR "Maximum coverage: $maxCov\n";
	print STDERR "Average coverage over $len reference-sequence bases: $avgCov\n";

	# Pass 2 draws the actual plot
	my $basesPerRow = min($opt->{maxBasesPerRow}, $len);
	print STDERR "Bases per row: $basesPerRow\n";

	ensureHeaderPrinted $g;

	print $outF <<THE_END;
% Things that we don't know until we've seen the data for a particular individual.
% Earlier-defined functions that use these names must not use bind.
/basesperrow $basesPerRow def
/maxcoverage $maxCov def
/w plotwidth basesperrow div def		% The width of a bar
/h plotheight maxcoverage div def		% The height of 1 coverage unit
THE_END
	
#	for (my $iRow = 0; $iRow < $nRows; ++$iRow) {
#		newRow $g;
#		print $outF "(", psString($inFn), ") title\n" if $iRow == 0;
#	}

	# Put all the logic for plotting the next point into a 2-arg closure to avoid code duplication.
	my $plotIt = sub {
		my ($i, $cov) = @_;

		if (($i - $opt->{startBase}) % $basesPerRow == 0) {		# Yes, we want this to trigger on the first base too.
			if ($i > $opt->{startBase}) {
				# Draw an average line on the previous row.
#				drawAverage $avgCov, $g;
				# If we get to here, then we know that we're on the 2nd or later row of
				# a multi-row plot, so the previous row must have had maxBasesPerRow bases.
				my $nBases = $opt->{maxBasesPerRow};
				drawRest $i, $nBases, $basesPerRow, $maxCov, $avgCov, $opt, $g;
			}

#			my $nBases = min($opt->{maxBasesPerRow}, $opt->{endBase} - $i + 1);
#			newRow $i, $nBases, $g;
			newRow $g;
			my $bn = $inFn;
			$bn =~ s/\.cov\z//i;			# Strip a ".cov" prefix if present
			print $outF "(", psString("$bn: " . sprintf "%d (%.2f%%) of %d bases covered", $nBasesCovered, $nBasesCovered / $len * 100, $len), ") title\n" if $i == $opt->{startBase};
		};

		if (($i - 1) % $opt->{xTicks} == 0) {		# Using "$i % $opt->{xTicks} == 1" doesn't work when xTicks = 1.
			print $outF "($i) xtick\n";		# Don't need to use psString() since $i is always directly representable.
		}

		print $outF "$cov p\n";
	};

	{
		local $_;
		open my $inF, "<", $inFn or die "Could not open '$inFn' for input: $!";
		my $oldLoc = $opt->{startBase} - 1;
		while (<$inF>) {
			my (undef, $loc, undef, $cov) = split /\t/;

			last if $loc > $opt->{endBase};
			if ($loc >= $opt->{startBase}) {
				# Every missing row means 0 coverage for that location.
				for (my $i = $oldLoc + 1; $i < $loc; ++$i) {
					$plotIt->($i, 0);
				}

				$plotIt->($loc, $cov);

				$oldLoc = $loc;
			}
		}
		close $inF;

		# Every missing row means 0 coverage for that location.
		for (my $i = $oldLoc + 1; $i <= $opt->{endBase}; ++$i) {
			$plotIt->($i, 0);
		}

#		drawAverage $avgCov, $g;
		my $nBases = $len % $basesPerRow;
		$nBases = $basesPerRow if !$nBases;
		drawRest $opt->{endBase} - $nBases + 1, $nBases, $basesPerRow, $maxCov, $avgCov, $opt, $g;
	}
}

# Main program
my $opt = { %$defaultOpt };			# Copies the elements
my @fNames;

my $syntax = <<THE_END;
Create PostScript coverage plots from coverage files.

Syntax:
  make_coverage_plots.pl [options] [cov_file1 [cov_file2 [cov_file3 ...]]] > covplot.ps

Options:
  -r range, --range range
     Specify the starting and ending nucleotide positions within the reference
	 sequence.  range must be of the form "x-y".  If not specified, the range
	 will be automatically determined from the minimum and maximum positions
	 mentioned in the coverage input files.
  
  -b bases, --max-bases-per-row bases
     Specify the maximum number of bases to fit on one row of a plot.
	 Default: $opt->{maxBasesPerRow}
  
  -t bases, --bases-between-ticks bases
     Specify the number of bases between tick-marks appearing on the x axis.
	 Default: $opt->{xTicks}
  
  -h, --help
     Show this help text.

  --
     Treat remaining command-line arguments as filenames, even if they look
     like options.

cov_file1 etc. should be text files that give the coverage level for each
nucleotide position within the reference sequence, in the format output by
samtools mpileup -D.  Nucleotide positions having zero coverage may be
omitted.  Lines must appear in increasing order of position.  Each file should
give the coverage information for a different individual.

The output file is in PostScript format, suitable for viewing with a program
like GhostView, or for printing directly to a PostScript printer.
THE_END

# Process command-line arguments
while (@ARGV) {
	$_ = shift;
	if ($_ eq '-r' || $_ eq '--range') {
		die "-r must be followed by a range.\n\n$syntax" if !@ARGV;
		$_ = shift;
		($opt->{startBase}, $opt->{endBase}) = /\A(\d+)-(\d+)\z/ or die "Range must be in format 'x-y'.\n\n$syntax";
	} elsif ($_ eq '-b' || $_ eq '--max-bases-per-row') {
		die "-b must be followed by the number of bases.\n\n$syntax" if !@ARGV;
		$opt->{maxBasesPerRow} = shift;
	} elsif ($_ eq '-t' || $_ eq '--bases-between-ticks') {
		die "-t must be followed by the number of bases between x-axis ticks.\n\n$syntax" if !@ARGV;
		$opt->{xTicks} = shift;
	} elsif ($_ eq '-h' || $_ eq '--help') {
		print $syntax;
		exit 0;
	} elsif ($_ eq '--') {
		push @fNames, @ARGV;
		last;
	} elsif (/\A-/) {
		die "Unrecognised option '$_'.\n\n$syntax";
	} else {
		push @fNames, $_;
	}
}

# The range can now be determined automatically if it is not specified.
#die "Must specify a range with -r.\n\n$syntax" if !defined $opt->{startBase} || !defined $opt->{endBase};

print STDERR "Maximum number of bases per plot row: $opt->{maxBasesPerRow}\n";

# A graphics context that writes PostScript to standard output.
my $g = {
	outF => \*STDOUT,
	maxRows => 4
};

foreach my $fn (@fNames) {
	plot $fn, $opt, $g;
}

finalise $g;
print STDERR "The end.\n";
