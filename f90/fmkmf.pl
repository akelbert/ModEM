#!/usr/bin/perl -w
# Copyright (c) The University of Edinburgh
# This is a utility to generate make files 
# for Fortran 90. It was originally in shell script and was re-written
# in perl for greater speed and (hopefully) portability.
# Initial tests suggest speed is 10x better than the sh version.
#
# Format modified by Anna Kelbert, 2005-2009.
#
# A basic makefile entry for bork.f90 would be
# bork.o:bork.f90
# <-tab->$(F90) -c bork.f90
#
# however if bork.f90 contains the line "use gunge" then 
# (A)
# the entry has to be 
# bork.o:bork.f90 garple.o <-- Forces bork to be recompiled if a module it 
# <-tab->$(F90) -c bork.f90                               uses is changed
# where garple.f90 is the program containing the line "module gunge
# (B)  
# The same type of entry has to be done for garple.f90
#
# We also need to generate an entry for the link step. If the main program
# was in baz.f90 then this should be 
# baz:baz.o bork.o.........
# <-tab->$(F90) -o baz baz.o bork.o .....
# The list of object files to be linked should have foo.o in it once 
# and only once for each foo.f90 that was compiled 

use File::Basename;

#-------------------------------------------------
# First check if the luser has any relevent environment vars set
#--------------------------------------------
if ( $ENV{FMKMF_F90} ) {
  #print "\# FMKMF_F90 set to $ENV{FMKMF_F90}\n";
  $f90=$ENV{FMKMF_F90};
} 
else {
  #print "\# FMKMF_F90 not set: using f90\n";
  $f90="f90";
}


if ( $ENV{FMKMF_SFTAG} ) {
  #print "\# FMKMF_SFTAG set to $ENV{FMKMF_SFTAG}\n";
  $sftag=$ENV{FMKMF_SFTAG};
} 
else {
  #print "\# FMKMF_SFTAG not set: using f90\n";
  $sftag="f90";
}

if ( $ENV{FMKMF_SPATH} ) {
  #print "\# FMKMF_SPATH set to $ENV{FMKMF_SPATH}\n";
  $spath=$ENV{FMKMF_SPATH};
} 
else {
  #print "\# FMKMF_SPATH not set: using . \n";
  $spath=".";
}

if ( $ENV{FMKMF_OPTIM} ) {
  #print "\# FMKMF_OPTIM set to $ENV{FMKMF_OPTIM}\n";
  $optim=$ENV{FMKMF_OPTIM};
} 
else {
  #print "\# FMKMF_OPTIM not set: using default optimization \n";
  $optim="-O2";
}

if ( $ENV{FMKMF_LINKOPTS} ) {
  #print "\# FMKMF_LINKOPTS set to $ENV{FMKMF_LINKOPTS}\n";
  $linkopts=$ENV{FMKMF_LINKOPTS};
} 
else {
  #print "\# FMKMF_LINKOPTS not set: using no link options \n";
  $linkopts=" ";
}

# By default, use the current directory for object files
$linkdir=".";

#------------------------------
# Done with environment variables. Now we need to process commandline args
# These supersede anything supplied via environment variables.
#------------------------------

my ($optiond)=0;

while (@ARGV){

  $arg=shift;
  if ($arg =~ /^-p$/){
    $spath=shift;
    #print "# Using search path $spath from cmd line\n";
  }
  if ($arg =~ /^-f90$/){
    $f90=shift;
    #print "# Using compile cmd $f90 from cmd line\n";
  }
  if ($arg =~ /^-tag$/){
    $sftag=shift;
    #print "# Using source file tag $sftag from cmd line\n";
  }
  if ($arg =~ /^-opt$/){
    $optim=shift;
    #print "# Using compiler optimization options $optim from cmd line\n";
  }
  if ($arg =~ /^-l$/){
    $linkopts=shift;
    #print "# Using Link options $linkopts from cmd line\n";
  }
  if ($arg =~ /^-d$/){
    $optiond=1
    #print "# Using debug option (full output on) from cmd line\n";
  }
  if ($arg =~ /^-o$/){
  	$linkdir=shift;
  	#print "# Using the default object file output directory\n";
  }

}

#-------------------------------------------
# Done processing command line args
#-------------------------------------------


@spath=split(/:/,$spath);


@global_outlines=();
@global_objlist=();
@global_modfiles=();

$mainprogfile=$arg;

# Generate a name for the executable file
$execfile=$mainprogfile;
$execfile=~s/\.${sftag}//;
$execfile=~s|.*/||;

# Output makefile header
print "# Makefile suited for building the $execfile program\n";
print "# Generated using: ./fmkmf.pl [OPTIONS] $mainprogfile > Makefile\n";
print "# with command line options\n";
if($optiond){
  print "# -d (extra debugging output)\n";
}
print "# -p $spath\n";
print "# -f90 $f90 (compiler)\n";
print "# -opt $optim (compiler optimisation)\n";
print "# -l $linkopts (linking options)\n";
print "# -o $linkdir (output directory for object files)\n\n";

print "#  Uncomment these lines to make program for Solaris OS\n";
print "# F90 = f90\n";
print "# FFLAGS = -dalign -g -C -w  -L/usr/local/lib\n";
print "# LIBS = -xlic_lib=sunperf\n";
print "#  Uncomment these lines to make program with g95\n";
print "# F90 = g95\n";
print "# FFLAGS = -O2\n";
print "# FFLAGS = -g -ftrace=frame -fbounds-check\n";
print "# MODULE = -fmod=\$(OBJDIR)\n";
print "# LIBS = -lblas -llapack\n";
print "#  Uncomment these lines to make program with Intel compiler\n";
print "# F90 = ifort\n";
print "# FFLAGS = -O3\n";
print "# FFLAGS = -g -debug all\n";
print "# MODULE = -module \$(OBJDIR)\n";
print "# LIBS = -lblas -llapack\n";
print "#  Uncomment these lines to make program with PGI compiler\n";
print "# F90 = pgf95\n";
print "# FFLAGS = -O3\n";
print "# FFLAGS = -g -Mprof=lines -Mbounds\n";
print "# MODULE = -module \$(OBJDIR)\n";
print "# LIBS = -L/usr/lib64 -lblas -llapack -lpgftnrtl -Mprof=lines\n";

if($optiond){
  print "# Main program is $mainprogfile \n" ;
}
# this subroutine (def below) does most of the work.
process_fsource($mainprogfile); 

# set some makefile . 

print "\n# ------------------Macro-Defs---------------------\n";

print "include Makefile.local\n";
print "OBJDIR = $linkdir\n";
print "F90 = $f90 \n";
print "FFLAGS = $optim\n";
print "MODULE = -module \$(OBJDIR)\n";
print "LIBS = $linkopts\n";


print "\n# -------------------End-macro-Defs---------------------------\n";

print "OBJ = @global_objlist \n\n";


print "\nall: $execfile \n";

# Generate makefile entry for the Link step
print "\n# Here is the link step \n";

print "$execfile: \$(OBJDIR) \$(OBJ) \n";
print "\t \$(F90) -o \$(OUTDIR)/$execfile \$(OBJ) \$(LIBS) \n";
# print "\trm -f *.mod \n";

print "\n# Here are the compile steps \n\n";

print "\$(OBJDIR): \n";
print "\tmkdir -p \$(OBJDIR)\n";

print STDOUT @global_outlines;

# Add an entry for make clean at the end of the make file.  this
# removes most of the garbage left around by most of the Fortran 90
# compilers I have tried.

print "\n# Type \" make clean \" to get rid of all object and module files \n";

print "clean: \n";
print "\tcd \$(OBJDIR); \\\n";
print "\trm -f *~ *.o *.obj *.mod *.d *.s00 *.dbg *.stackdump \\\n";
print "\t`find . -mindepth 1 -name \"*~\"` \n\n";
print "cleanall: clean \n";
print "\trm -f \$(OUTDIR)/$execfile \n\n";
print "src: clean \n";
print "\ttar cvfz \$(ARCHIVE).tgz * \n";
print "\n  \n";

# End of main program 

##############################################
# Here is the subroutine that generates the compile entries in the makefile
# These end up in the global array @global_outlines. The magic part is 
# that this subroutine calls itself recursively.
##############################################
sub process_fsource {
  
  my $mainprogfile=$_[0];
  if($optiond){
    print"# process_fsource called with arg $mainprogfile \n";
  }
  open( MAINPROG, $mainprogfile) or 
    die "Can't find main program file $mainprogfile: $! \n";
  
  # Read through Fortran source looking for USE statements
  # There should be nothing but whitespace before the USE. Sloppily, 
  # we allow tabs, although the standard (IIRC) does not
  my @modulelist=();
  my @includelist=();
  while ($line=<MAINPROG>) { 
    if ($line =~ /^[ \t]*use (\w+)/i ) { # line matches regexp between / /
      if($optiond){
		print "# $mainprogfile Uses Module $1\n";
      }
      @modulelist=(@modulelist,$1);
    } elsif  ($line =~ /^[ \t]*include "(\S+)"/i ){
      if($optiond){
		print "# $mainprogfile Includes $1\n";
      }
      my $includefile = $1;
      my $mainprogpath = dirname($mainprogfile); 
      @includelist=(@includelist,"$mainprogpath/$includefile");    	
    }
  }
  
  close(MAINPROG);

  if($optiond){
    print "# Full list of modules in $mainprogfile: @modulelist \n";
    if (@includelist > 0) {
    	print "# Full list of includes in $mainprogfile: @includelist \n";
    }
  }
  # Find which file each module is in.
  
  
  
 my @modfiles=();
 MODLOOP:foreach $module (@modulelist){
    foreach $directory (@spath){
      # print "# Looking in directory $directory\n";
      opendir( DIRHANDLE, $directory) or die 
	"Can't open directory $directory : $! \n";
      @sourcefiles=grep /\.${sftag}\Z/, sort(readdir(DIRHANDLE));
    foreach $sourcefile (@sourcefiles){
      $pathsourcefile="$directory/$sourcefile";
      #print "\# Checking $pathsourcefile\n";
      open( SOURCEFILE, "$pathsourcefile") or 
	die "Can't find source file $pathsourcefile: $! \n";
      while ($line=<SOURCEFILE>){
	if ($line =~ /^ *module (\w+)/i ){
	  if($1 =~ /^$module$/i){
	    if($optiond){
	      print "# Uses $module which is in $pathsourcefile\n";
	    }
	    @modfiles=(@modfiles,$pathsourcefile);
	    
	    if (grep (/$pathsourcefile/,@global_modfiles )){
	      if($optiond){
		print "# $pathsourcefile already in list\n";
	      }
	    }
	    else {
	      @global_modfiles=(@global_modfiles,$pathsourcefile);
	      process_fsource($pathsourcefile);

	    }
	    # We found this module -- go on to the next one
	    close (SOURCEFILE);
	    next MODLOOP;	    
	  }
	}
      }
      close( SOURCEFILE );
    }
  }
  # exhausted source files
  print STDERR "Couldn't find source file for module $module\n";
}

# name of file we want to make
$objfile=$mainprogfile;
# replace source file name with .o
$objfile=~s/\.${sftag}/\.o/;
# strip path so object files go in current dir
$objfile=~s|.*/||;
# now add the user-defined path to the object files
$objfile="\$(OBJDIR)/$objfile";
@global_objlist=(@global_objlist,$objfile);
# list of dependencies
@objlist=();
foreach  $mf (@modfiles) { 
  $obj=$mf;
  # replace source file name with .o
  $obj=~s/\.${sftag}/\.o/;
  # strip path so object files go in current dir
  $obj=~s|.*/||;
  # now add the user-defined path to the object file
  $obj="\$(OBJDIR)/$obj";
  @objlist=(@objlist,$obj);
}

@global_outlines=(@global_outlines,"\n$objfile:$mainprogfile @objlist @includelist\n");
@global_outlines=(@global_outlines,"\t \$(F90) -c \$(MODULE) \$(FFLAGS) $mainprogfile -o $objfile\n");

#if (@includelist > 0) {
#	@global_outlines=(@global_outlines,"\n$mainprogfile: @includelist \n");
#}

}
