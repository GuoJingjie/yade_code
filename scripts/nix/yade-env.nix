with (import <nixpkgs> {});

( let
minieigen = pkgs.python37Packages.buildPythonPackage rec{
      name = "minieigen";

      src = pkgs.fetchFromGitHub {
        owner = "eudoxos";
        repo = "minieigen";
        rev = "1e992b1452638e636b6681a9ab17f6bbb51a727d";
        sha256 = "1a0c01kzgiy3p00ifids0qwqdkyaly5iiy5mf76ymn0hzcwldnfx";
      };

      buildInputs = [ unzip python37Packages.boost eigen ];

      patchPhase = ''
        sed -i "s/^.*libraries=libraries.//g" setup.py 
      '';

      preConfigure = ''
        export LDFLAGS="-L${eigen.out}/lib -L${python37Packages.boost.out} -lboost_python37"
        export CFLAGS="-I${eigen.out}/include/eigen3"
      '';
};

in 

{ yade-env = pkgs.python37.buildEnv.override rec{

        extraLibs = with pkgs.python37Packages;[
                        pygments mpi4py pexpect decorator numpy xlib
                        ipython ipython_genutils traitlets pygraphviz
                        six minieigen ipython future matplotlib tkinter pillow pkgs.cmake
                      ] ;
        ignoreCollisions = true;
    };
}
)
