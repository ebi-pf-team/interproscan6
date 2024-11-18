class Output {
    static formatCathAccession (String modelAccession) {
        return "G3DSA:" + modelAccession.replace("-", ":")
    }

    static convertDbName (String memberDB) {
        switch (memberDB) {
            case "antifam":
                return "AntiFam"
            case "cdd":
                return "CDD"
            case "coils":
                return "Coils"
            case "cathfunfam":
                return "FunFam"
            case "cathgene3d":
                return "Gene3D"
            case "hamap":
                return "Hamap"
            case "mobidblite":
                return "MobiDBLite"
            case "ncbifam":
                return "NCBIfam"
            case "panther":
                return "PANTHER"
            case "pfam":
                return "Pfam"
            case "pirsf":
                return "PIRSF"
            case "pirsr":
                return "PIRSR"
            case "phobius":
                return "Phobius"
            case "prints":
                return "PRINTS"
            case "prositepatterns":
                return "ProSitePatterns"
            case "prositeprofiles":
                return "ProSiteProfiles"
            case "sfld":
                return "SFLD"
            case "signalp":
                return "SignalP"
            case "signalp_euk":
                return "SignalP_EUK"
            case "smart":
                return "SMART"
            case "superfamily":
                return "SUPERFAMILY"
            case "tmhmm":
                return "DeepTMHMM"
        }
    }
}
