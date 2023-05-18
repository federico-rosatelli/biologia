package structures

type Taxonomy struct {
	TaxID          string `bson:"TaxId"`
	ScientificName string `bson:"ScientificName"`
	ParentTaxID    string `bson:"ParentTaxId"`
	Rank           string `bson:"Rank"`
	Division       string `bson:"Division"`
	GeneticCode    struct {
		GCID   string `bson:"GCId"`
		GCName string `bson:"GCName"`
	} `bson:"GeneticCode"`
	MitoGeneticCode struct {
		MGCID   string `bson:"MGCId"`
		MGCName string `bson:"MGCName"`
	} `bson:"MitoGeneticCode"`
	Lineage   string `bson:"Lineage"`
	LineageEx []struct {
		TaxID          string `bson:"TaxId"`
		ScientificName string `bson:"ScientificName"`
		Rank           string `bson:"Rank"`
	} `bson:"LineageEx"`
	CreateDate string `bson:"CreateDate"`
	UpdateDate string `bson:"UpdateDate"`
	PubDate    string `bson:"PubDate"`
}

type TaxonomyTree struct {
	TaxId          string `bson:"TaxId"`
	Rank           string `bson:"Rank"`
	ScientificName string `bson:"ScientificName"`
	SubClasses     []struct {
		TaxID          string `bson:"TaxId"`
		ScientificName string `bson:"ScientificName"`
		Rank           string `bson:"Rank"`
	} `bson:"SubClasses"`
}

type OrganismTable struct {
	ScientificName string
	TaxId          string
	QtyNucleotides int
	QtyProteins    int
	QtyProducts    int
	Genomes        []string
	Annotations    []string
	Trascriptome   []string
	SraWgs         int
	SraTran        int
}

type Nucleotide struct {
	GBSeqAccessionVersion string `bson:"GBSeq_accession-version"`
	GBSeqComment          string `bson:"GBSeq_comment"`
	GBSeqCreateDate       string `bson:"GBSeq_create-date"`
	GBSeqDefinition       string `bson:"GBSeq_definition"`
	GBSeqDivision         string `bson:"GBSeq_division"`
	GBSeqFeatureTable     []struct {
		GBFeatureIntervals []struct {
			GBIntervalAccession string `bson:"GBInterval_accession"`
			GBIntervalFrom      string `bson:"GBInterval_from"`
			GBIntervalTo        string `bson:"GBInterval_to"`
		} `bson:"GBFeature_intervals"`
		GBFeatureKey      string `bson:"GBFeature_key"`
		GBFeatureLocation string `bson:"GBFeature_location"`
		GBFeatureQuals    []struct {
			GBQualifierName  string `bson:"GBQualifier_name"`
			GBQualifierValue string `bson:"GBQualifier_value"`
		} `bson:"GBFeature_quals"`
		GBFeaturePartial3 string `bson:"GBFeature_partial3,omitempty"`
		GBFeaturePartial5 string `bson:"GBFeature_partial5,omitempty"`
	} `bson:"GBSeq_feature-table"`
	GBSeqLength           string   `bson:"GBSeq_length"`
	GBSeqLocus            string   `bson:"GBSeq_locus"`
	GBSeqMoltype          string   `bson:"GBSeq_moltype"`
	GBSeqOrganism         string   `bson:"GBSeq_organism"`
	GBSeqOtherSeqids      []string `bson:"GBSeq_other-seqids"`
	GBSeqPrimaryAccession string   `bson:"GBSeq_primary-accession"`
	GBSeqReferences       []struct {
		GBReferenceAuthors   []string `bson:"GBReference_authors"`
		GBReferenceJournal   string   `bson:"GBReference_journal"`
		GBReferencePosition  string   `bson:"GBReference_position"`
		GBReferenceReference string   `bson:"GBReference_reference"`
		GBReferenceTitle     string   `bson:"GBReference_title"`
	} `bson:"GBSeq_references"`
	GBSeqSource       string `bson:"GBSeq_source"`
	GBSeqStrandedness string `bson:"GBSeq_strandedness"`
	GBSeqTaxonomy     string `bson:"GBSeq_taxonomy"`
	GBSeqTopology     string `bson:"GBSeq_topology"`
	GBSeqUpdateDate   string `bson:"GBSeq_update-date"`
}

type NucleotideBasic struct {
	ScientificName string `bson:"ScientificName"`
	GBSeq_locus    []struct {
		GBSeq_locus string `bson:"GBSeq_locus"`
	} `bson:"GBSeq_locus"`
}
