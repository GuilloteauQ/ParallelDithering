let
    pkgs = import <nixpkgs> {};
    stdenv = pkgs.stdenv;
in with pkgs; {
    presentation = stdenv.mkDerivation {
        name = "report";
        buildInputs = [
            R
            texlive.combined.scheme-full
            pandoc
            rPackages.knitr
            rPackages.rmarkdown
            rPackages.ggplot2
            rPackages.dplyr
        ];
        src = ./.;
        installPhase = ''
          mkdir $out
          Rscript -e "rmarkdown::render('main.Rmd', 'pdf_document')"
          cp main.pdf $out
        '';
    };
}
