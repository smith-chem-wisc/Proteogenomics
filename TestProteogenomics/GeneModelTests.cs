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
        public void Setup()
        {
            genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));

            string proteinFasta = "Homo_sapiens.GRCh38.pep.all.fa";
            string chr5Fasta = "Homo_sapiens.GRCh38.dna.chromosome.5.fa";
            string chr7Fasta = "Homo_sapiens.GRCh38.dna.chromosome.7.fa";
            string chr14Fasta = "Homo_sapiens.GRCh38.dna.chromosome.14.fa";
            string chr19Fasta = "Homo_sapiens.GRCh38.dna.chromosome.19.fa";
            string chr20Fasta = "Homo_sapiens.GRCh37.73.dna.chromosome.20.fa";
            string chr21Fasta = "Homo_sapiens.GRCh37.73.dna.chromosome.21.fa";
            string chr22Fasta = "Homo_sapiens.GRCh37.73.dna.chromosome.22.fa";
            string chrMTFasta = "Homo_sapiens.GRCh38.dna.chromosome.MT.fa";
            string[] gunzippedFiles = new[] 
            {
                proteinFasta,
                chr5Fasta,
                chr7Fasta,
                chr14Fasta,
                chr19Fasta,
                chr20Fasta,
                chr21Fasta,
                chr22Fasta,
                chrMTFasta,
            };

            // download files
            if (!gunzippedFiles.All(f => File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, f)) && new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, f)).Length > 0))
            {
                using (WebClient Client = new WebClient())
                {
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, proteinFasta + ".gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, chr5Fasta + ".gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, chr7Fasta + ".gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, chr14Fasta + ".gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, chr19Fasta + ".gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.73.dna.chromosome.20.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, chr20Fasta + ".gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.73.dna.chromosome.21.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, chr21Fasta + ".gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.73.dna.chromosome.22.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, chr22Fasta + ".gz"));
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz", Path.Combine(TestContext.CurrentContext.TestDirectory, chrMTFasta + ".gz"));
                }
            }

            // decompress them
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

            // combine 20, 21, and 22 for a test
            string fasta202122Path = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            if (!File.Exists(fasta202122Path))
            {
                using (var f = File.Create(fasta202122Path))
                {
                    foreach (var fasta in new[] { chr20Fasta, chr21Fasta, chr22Fasta })
                    {
                        using (var stream = new FileStream(Path.Combine(TestContext.CurrentContext.TestDirectory, fasta), FileMode.Open))
                        {
                            stream.CopyTo(f);
                        }
                    }
                }
            }
        }

        [Test]
        public void GtfBasics()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.AreEqual(165, geneModel.Genes.SelectMany(g => g.Transcripts).Count());
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Translate(true)).ToList();
        }

        [Test]
        public void GffBasics()
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
        public void GffAppliedToOther()
        {
            string referenceGff = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3");
            string alternateGff = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_pacbio.gtf");
            GeneModel r = new GeneModel(genome, referenceGff);
            GeneModel a = new GeneModel(genome, alternateGff);
            a.CreateCDSFromAnnotatedStartCodons(r);
            List<Protein> proteins = a.Genes.SelectMany(g => g.Translate(true)).ToList();

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
        public void OutputGtfFromGeneModel()
        {
            string referenceGff = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3");
            string alternateGff = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_pacbio.gtf");
            GeneModel r = new GeneModel(genome, referenceGff);
            GeneModel a = new GeneModel(genome, alternateGff);
            a.CreateCDSFromAnnotatedStartCodons(r);
            a.PrintToGTF(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_pacbio_merged.gtf"));
        }
    }
}