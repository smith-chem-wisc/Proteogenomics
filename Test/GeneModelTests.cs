using Bio;
using Bio.Extensions;
using Bio.IO.FastA;
using NUnit.Framework;
using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Net;

namespace TestProteogenomics
{
    [TestFixture]
    public class GeneModelTests
    {
        private Genome genome;

        [OneTimeSetUp]
        public void setup()
        {
            genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));

            string proteinFasta = "Homo_sapiens.GRCh38.pep.all.fa";
            string chr5Fasta = "Homo_sapiens.GRCh38.dna.chromosome.5.fa";
            string chr7Fasta = "Homo_sapiens.GRCh38.dna.chromosome.7.fa";
            string chr14Fasta = "Homo_sapiens.GRCh38.dna.chromosome.14.fa";
            string chr19Fasta = "Homo_sapiens.GRCh38.dna.chromosome.19.fa";
            string chrMTFasta = "Homo_sapiens.GRCh38.dna.chromosome.MT.fa";
            string[] gunzippedFiles = new[] { proteinFasta, chr5Fasta, chr7Fasta, chr14Fasta, chr19Fasta, chrMTFasta };
            if (!gunzippedFiles.All(f => File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, f)) && new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, f)).Length > 0))
            {
                using (WebClient Client = new WebClient())
                {
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.pep.all.fa.gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, "Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"));
                }
            }

            foreach (var gunzippedFile in gunzippedFiles)
            {
                if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, gunzippedFile)) || new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, gunzippedFile)).Length == 0)
                {
                    using (FileStream stream = new FileStream(Path.Combine(TestContext.CurrentContext.TestDirectory, gunzippedFile + ".gz"), FileMode.Open))
                    using (GZipStream gunzip = new GZipStream(stream, CompressionMode.Decompress))
                    using (var f = File.Create(Path.Combine(TestContext.CurrentContext.TestDirectory, gunzippedFile)))
                    {
                        gunzip.CopyTo(f);
                    }
                }
            }
        }

        [Test]
        public void CountReads()
        {
            Assert.AreEqual(3970, new FastqProperties(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper.fastq")).ReadCount);
        }

        [Test]
        public void gtfBasics()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.AreEqual(165, geneModel.Genes.SelectMany(g => g.Transcripts).Count());
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Translate(true)).ToList();
        }

        [Test]
        public void gffBasics()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.AreEqual(148, geneModel.Genes.SelectMany(g => g.Transcripts).Count());
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Translate(true)).ToList();

            //Forward strand, single coding region
            Assert.AreEqual("ENSP00000334393", proteins[0].Accession);
            Assert.AreEqual(
                "MVTEFIFLGLSDSQELQTFLFMLFFVFYGGIVFGNLLIVITVVSDSHLHSPMYFLLANLSLIDLSLSSVTAPKMITDFFSQRKVISFKGCLVQIFLLHFFGGSEMVILIAMGFDRYIAICKPLHYTTIMCGNACVGIMAVTWGIGFLHSVSQLAFAVHLLFCGPNEVDSFYCDLPRVIKLACTDTYRLDIMVIANSGVLTVCSFVLLIISYTIILMTIQHRPLDKSSKALSTLTAHITVVLLFFGPCVFIYAWPFPIKSLDKFLAVFYSVITPLLNPIIYTLRNKDMKTAIRQLRKWDAHSSVKF",
                proteins[0].BaseSequence);

            //Reverse strand, single coding region
            Assert.AreEqual("ENSP00000473460", proteins[1].Accession);
            Assert.AreEqual(
                "TSLWTPQAKLPTFQQLLHTQLLPPSGLFRPSSCFTRAFPGPTFVSWQPSLARFLPVSQQP" +
                "RQAQVLPHTGLSTSSLCLTVASPRPTPVPGHHLRAQNLLKSDSLVPTAASWWPMKAQNLL" +
                "KLTCPGPAPASCQRLQAQPLPHGGFSRPTSSSWLGLQAQLLPHNSLFWPSSCPANGGQCR" +
                "PKTSSSQTLQAHLLLPGGINRPSFDLRTASAGPALASQGLFPGPALASWQLPQAKFLPAC" +
                "QQPQQAQLLPHSGPFRPNL",
                proteins[1].BaseSequence);
        }

        [Test]
        public void gffAppliedToOther()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            GeneModel additional = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_pacbio.gff3"));
            geneModel.Merge(additional);
            List<Protein> proteins = additional.Genes.SelectMany(g => g.Translate(true)).ToList();

            //Forward strand, single coding region
            Assert.AreEqual("PB2015.1.1", proteins[0].Accession);
            Assert.AreEqual(
                "MVTEFIFLGLSDSQELQTFLFMLFFVFYGGIVFGNLLIVITVVSDSHLHSPMYFLLANLSLIDLSLSSVTAPKMITDFFSQRKVISFKGCLVQIFLLHFFGGSEMVILIAMGFDRYIAICKPLHYTTIMCGNACVGIMAVTWGIGFLHSVSQLAFAVHLLFCGPNEVDSFYCDLPRVIKLACTDTYRLDIMVIANSGVLTVCSFVLLIISYTIILMTIQHRPLDKSSKALSTLTAHITVVLLFFGPCVFIYAWPFPIKSLDKFLAVFYSVITPLLNPIIYTLRNKDMKTAIRQLRKWDAHSSVKF",
                proteins[0].BaseSequence);

            //Reverse strand, single coding region
            Assert.AreEqual("PB2015.2.1", proteins[1].Accession);
            Assert.AreEqual(
                "TSLWTPQAKLPTFQQLLHTQLLPPSGLFRPSSCFTRAFPGPTFVSWQPSLARFLPVSQQP" +
                "RQAQVLPHTGLSTSSLCLTVASPRPTPVPGHHLRAQNLLKSDSLVPTAASWWPMKAQNLL" +
                "KLTCPGPAPASCQRLQAQPLPHGGFSRPTSSSWLGLQAQLLPHNSLFWPSSCPANGGQCR" +
                "PKTSSSQTLQAHLLLPGGINRPSFDLRTASAGPALASQGLFPGPALASWQLPQAKFLPAC" +
                "QQPQQAQLLPHSGPFRPNL",
                proteins[1].BaseSequence);
        }

        [Test]
        public void KaryotypicOrder()
        {
            Genome headers = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "headers.fa"));
            var seqs = headers.KaryotypicOrder();
            Assert.IsTrue(seqs[0].FriendlyName == "chr1" && seqs[1].FriendlyName == "chr2");
            Assert.IsTrue(seqs[22].FriendlyName == "chrX" && seqs[23].FriendlyName == "chrY" && seqs[24].FriendlyName == "chrM");
        }

        [Test]
        public void KaryotypicOrderShort()
        {
            Genome headers = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "headersShort.fa"));
            var seqs = headers.KaryotypicOrder();
            Assert.IsTrue(seqs[0].FriendlyName == "chr9" && seqs[1].FriendlyName == "chr20");
        }

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

            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory,  "Homo_sapiens.GRCh38.dna.chromosome.7.fa"));
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
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory,  "Homo_sapiens.GRCh38.dna.chromosome.5.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr5_selenocysteineContaining.gff3"));
            ProteinAnnotation.GetImportantProteinAccessions(Path.Combine(TestContext.CurrentContext.TestDirectory,  "Homo_sapiens.GRCh38.pep.all.fa"),
                out Dictionary<string, string> p, out HashSet<string> bad, out Dictionary<string, string> se);
            List<Protein> proteins = geneModel.Translate(true, bad, se).ToList();
            Assert.AreEqual("MWRSLGLALALCLLPSGGTESQDQSSLCKQPPAWSIRDQDPMLNSNGSVTVVALLQASUYLCILQASKLEDLRVKLKKEGYSNISYIVVNHQGISSRLKYTHLKNKVSEHIPVYQQEENQTDVWTLLNGSKDDFLIYDRCGRLVYHLGLPFSFLTFPYVEEAIKIAYCEKKCGNCSLTTLKDEDFCKRVSLATVDKTVETPSPHYHHEHHHNHGHQHLGSSELSENQQPGAPNAPTHPAPPGLHHHHKHKGQHRQGHPENRDMPASEDLQDLQKKLCRKRCINQLLCKLPTDSELAPRSUCCHCRHLIFEKTGSAITUQCKENLPSLCSUQGLRAEENITESCQURLPPAAUQISQQLIPTEASASURUKNQAKKUEUPSN",
                proteins[0].BaseSequence);
        }

        [Test]
        public void TranslateMTSeq()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory,  "Homo_sapiens.GRCh38.dna.chromosome.MT.fa"));

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