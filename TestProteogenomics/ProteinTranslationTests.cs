using Bio;
using Bio.Extensions;
using Bio.IO.FastA;
using NUnit.Framework;
using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TestProteogenomics
{
    [TestFixture]
    public class ProteinTranslationTests
    {
        [Test]
        public void TranslateReverseStrand()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_one_transcript_reverse.gtf"));
            List<Protein> proteins_wo_variant = geneModel.Translate(true).ToList();
            Assert.AreEqual("FFYFIIWSLTLLPRAGLELLTSSDPPASASQSVGITGVSHHAQ",
                proteins_wo_variant[0].BaseSequence);
        }

        [Test]
        public void TranslateAnotherReverseStrand()
        {
            // See http://useast.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000233306;r=7:38362864-38363518;t=ENST00000426402

            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.7.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr7_one_transcript_reverse.gtf"));
            List<Protein> proteins = geneModel.Translate(true).ToList();
            Assert.AreEqual("MQWALAVLLAFLSPASQKSSNLEGRTKSVIRQTGSSAEITCDLAEGSNGYIHWYLHQEGKAPQRLQYYDSYNSKVVLESGVSPGKYYTYASTRNNLRLILRNLIENDFGVYYCATWDG",
                proteins[0].BaseSequence);
        }

        [Test]
        public void TranslateHardReverseStrand()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.14.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "HardReverseStrand", "reverse.gff3"));
            List<Protein> proteins = geneModel.Translate(true).ToList();
            ISequence codingSequence = new FastAParser().Parse(Path.Combine(TestContext.CurrentContext.TestDirectory, "HardReverseStrand", "codingSeq.fa")).First();
            Assert.AreEqual(SequenceExtensions.ConvertToString(codingSequence),
                SequenceExtensions.ConvertToString(geneModel.Genes[0].Transcripts[0].RetrieveCodingSequence()));
            Assert.AreEqual("MNLQAQPKAQNKRKRCLFGGQEPAPKEQPPPLQPPQQSIRVKEEQYLGHEGPGGAVSTSQ" +
                "PVELPPPSSLALLNSVVYGPERTSAAMLSQQVASVKWPNSVMAPGRGPERGGGGGVSDSS" +
                "WQQQPGQPPPHSTWNCHSLSLYSATKGSPHPGVGVPTYYNHPEALKREKAGGPQLDRYVR" +
                "PMMPQKVQLEVGRPQAPLNSFHAAKKPPNQSLPLQPFQLAFGHQVNRQVFRQGPPPPNPV" +
                "AAFPPQKQQQQQQPQQQQQQQQAALPQMPLFENFYSMPQQPSQQPQDFGLQPAGPLGQSH" +
                "LAHHSMAPYPFPPNPDMNPELRKALLQDSAPQPALPQVQIPFPRRSRRLSKEGILPPSAL" +
                "DGAGTQPGQEATGNLFLHHWPLQQPPPGSLGQPHPEALGFPLELRESQLLPDGERLAPNG" +
                "REREAPAMGSEEGMRAVSTGDCGQVLRGGVIQSTRRRRRASQEANLLTLAQKAVELASLQ" +
                "NAKDGSGSEEKRKSVLASTTKCGVEFSEPSLATKRAREDSGMVPLIIPVSVPVRTVDPTE" +
                "AAQAGGLDEDGKGPEQNPAEHKPSVIVTRRRSTRIPGTDAQAQAEDMNVKLEGEPSVRKP" +
                "KQRPRPEPLIIPTKAGTFIAPPVYSNITPYQSHLRSPVRLADHPSERSFELPPYTPPPIL" +
                "SPVREGSGLYFNAIISTSTIPAPPPITPKSAHRTLLRTNSAEVTPPVLSVMGEATPVSIE" +
                "PRINVGSRFQAEIPLMRDRALAAADPHKADLVWQPWEDLESSREKQRQVEDLLTAACSSI" +
                "FPGAGTNQELALHCLHESRGDILETLNKLLLKKPLRPHNHPLATYHYTGSDQWKMAERKL" +
                "FNKGIAIYKKDFFLVQKLIQTKTVAQCVEFYYTYKKQVKIGRNGTLTFGDVDTSDEKSAQ" +
                "EEVEVDIKTSQKFPRVPLPRRESPSEERLEPKREVKEPRKEGEEEVPEIQEKEEQEEGRE" +
                "RSRRAAAVKATQTLQANESASDILILRSHESNAPGSAGGQASEKPREGTGKSRRALPFSE" +
                "KKKKTETFSKTQNQENTFPCKKCGR",
                proteins[0].BaseSequence);
        }

        [Test]
        public void TranslateSelenocysteineContaining()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.5.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr5_selenocysteineContaining.gff3"));
            ProteinAnnotation.GetImportantProteinAccessions(Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.pep.all.fa"),
                out Dictionary<string, string> p, out HashSet<string> bad, out Dictionary<string, string> se);
            List<Protein> proteins = geneModel.Translate(true, bad, se).ToList();
            Assert.AreEqual("MWRSLGLALALCLLPSGGTESQDQSSLCKQPPAWSIRDQDPMLNSNGSVTVVALLQASUYLCILQASKLEDLRVKLKKEGYSNISYIVVNHQGISSRLKYTHLKNKVSEHIPVYQQEENQTDVWTLLNGSKDDFLIYDRCGRLVYHLGLPFSFLTFPYVEEAIKIAYCEKKCGNCSLTTLKDEDFCKRVSLATVDKTVETPSPHYHHEHHHNHGHQHLGSSELSENQQPGAPNAPTHPAPPGLHHHHKHKGQHRQGHPENRDMPASEDLQDLQKKLCRKRCINQLLCKLPTDSELAPRSUCCHCRHLIFEKTGSAITUQCKENLPSLCSUQGLRAEENITESCQURLPPAAUQISQQLIPTEASASURUKNQAKKUEUPSN",
                proteins[0].BaseSequence);
        }

        [Test]
        public void TranslateMTSeq()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.MT.fa"));

            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chrM_one_transcript_reverse.gtf"));
            List<Protein> proteins = geneModel.Translate(true).ToList();
            Assert.AreEqual("MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYGLLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGLLFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYE" +
                "VTLAIILLSTLLMSGSFNLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTAYPRFRYDQLMHLLWKNFLPLTLALLMWYVSMPITISSIPPQT",
                proteins[0].BaseSequence);

            geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chrM_one_transcript_reverse2.gtf"));
            proteins = geneModel.Translate(true).ToList();
            Assert.AreEqual("MNPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVLTKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAMAMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPS" +
                "LNVSLLLTLSILSIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLLLNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSLIIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKM" +
                "KWQFEHTKPTPFLPTLIALTTLLLPISPFMLMIL",
                proteins[0].BaseSequence);
        }
    }
}