with import <nixpkgs>{};

let
  version = "1.0";
in {
  dithering = stdenv.mkDerivation rec {
    name = "dithering-${version}";
    src = ./.;
    buildInputs = [
      gcc
      coreutils
    ];
    installPhase = ''
      make parallel_dithering
      cp parallel_dithering ./bin
      cp -r ./bin $out
    '';
  };
}
