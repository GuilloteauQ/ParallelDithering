with import <nixpkgs>{};

let
  version = "1.0";
in {
  mpi_dithering = stdenv.mkDerivation rec {
    name = "mpi_dithering-${version}";
    src = ./.;
    buildInputs = [
      gcc
      coreutils
      openmpi
    ];
    installPhase = ''
      make dithering_bw
      cp dithering_bw ./bin
      cp -r ./bin $out
    '';
  };
}
