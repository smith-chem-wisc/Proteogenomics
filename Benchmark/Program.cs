using Proteogenomics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Net;
using UsefulProteomicsDatabases;

namespace Benchmark
{
    internal class Program
    {
        public static void Main(string[] args)
        {
            PacBioCds();
            //SameProteins();
        }

        /// <summary>
        /// Annotates PacBio transcript model for MCF7 with start codons in the reference model
        /// </summary>
        private static void PacBioCds()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            Genome genome = new Genome(@"E:\ProjectsActive\MCF7PacBio\Homo_sapiens.GRCh37.73.dna.primary_assembly.fa");
            string referenceGff = @"E:\ProjectsActive\MCF7PacBio\Homo_sapiens.GRCh37.73.gtf";
            string alternateGff = @"E:\ProjectsActive\MCF7PacBio\IsoSeq_MCF72015edition_polished.unimapped.ensembl.unimapped.gff";
            GeneModel r = new GeneModel(genome, referenceGff);
            GeneModel a = new GeneModel(genome, alternateGff);
            a.CreateCDSFromAnnotatedStartCodons(r);
            a.PrintToGTF(@"E:\ProjectsActive\MCF7PacBio\CDSAnnotated_IsoSeq_MCF7_2015edition_polished.unimapped.gff");

            stopwatch.Stop();
            Console.WriteLine("Finished checking that all proteins are the same.");
            Console.WriteLine("Time elapsed: " + stopwatch.Elapsed.Minutes.ToString() + " minutes and " + stopwatch.Elapsed.Seconds.ToString() + " seconds.");
            Console.WriteLine("Result: there are " + a.Genes.Sum(g => g.Transcripts.Count) + " PacBio transcript isoforms");
            Console.WriteLine("Result: " + a.Genes.Sum(g => g.Transcripts.Count(t => t.IsProteinCoding())) + " PacBio transcript isoforms are new annotated as protein coding");
            Console.WriteLine("Press any key to continue...");
            Console.ReadKey();
        }

        /// <summary>
        /// Times and checks that all proteins in the pep.all.fasta protein fasta file are the same as are output by this library
        /// </summary>
        private static void SameProteins()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            // download and decompress references
            string genomeFasta = "Homo_sapiens.GRCh38.dna.primary_assembly.fa";
            string geneModelFile = "Homo_sapiens.GRCh38.81.gff3";
            string proteinFasta = "Homo_sapiens.GRCh38.pep.all.fa";
            string[] gunzippedFiles = new[] { genomeFasta, geneModelFile, proteinFasta };
            if (!gunzippedFiles.All(f => File.Exists(f) && new FileInfo(f).Length > 0))
            {
                using (WebClient Client = new WebClient())
                {
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz");
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz", "Homo_sapiens.GRCh38.81.gff3.gz");
                    Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz", "Homo_sapiens.GRCh38.pep.all.fa.gz");
                }
            }

            foreach (var gunzippedFile in gunzippedFiles)
            {
                if (!File.Exists(gunzippedFile) || new FileInfo(gunzippedFile).Length == 0)
                {
                    using (FileStream stream = new FileStream(gunzippedFile + ".gz", FileMode.Open))
                    using (GZipStream gunzip = new GZipStream(stream, CompressionMode.Decompress))
                    using (var f = File.Create(gunzippedFile))
                    {
                        gunzip.CopyTo(f);
                    }
                }
            }

            ProteinAnnotation.GetImportantProteinAccessions(proteinFasta, out Dictionary<string, string> proteinAccessionSequence, out HashSet<string> bad, out Dictionary<string, string> se);

            Genome genome = new Genome(genomeFasta);
            GeneModel geneModel = new GeneModel(genome, geneModelFile);
            List<Protein> geneBasedProteins = geneModel.Translate(true, bad, se);
            List<Protein> pepDotAll = ProteinDbLoader.LoadProteinFasta(proteinFasta, true, DecoyType.None, false,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblGeneNameRegex, null, out List<string> errors);
            Dictionary<string, string> accSeq = geneBasedProteins.ToDictionary(p => p.Accession, p => p.BaseSequence);
            stopwatch.Stop();

            bool allAreEqual = true;
            foreach (Protein p in pepDotAll)
            {
                // now handled with the badAccessions // && !p.BaseSequence.Contains('*') && !seq.Contains('*') && !p.BaseSequence.Contains('X'))
                if (accSeq.TryGetValue(p.Accession, out string seq))
                {
                    if (p.BaseSequence != seq)
                    {
                        allAreEqual = false;
                        break;
                    }
                }
            }

            stopwatch.Stop();
            Console.WriteLine("Finished checking that all proteins are the same.");
            Console.WriteLine("Time elapsed: "+ stopwatch.Elapsed.Minutes.ToString() + " minutes and " + stopwatch.Elapsed.Seconds.ToString() + " seconds.");
            Console.WriteLine("Result: all proteins are " + (allAreEqual ? "" : "not ") + "equal ");
            Console.WriteLine("Press any key to continue...");
            Console.ReadKey();

            foreach (var file in Directory.GetFiles(Environment.CurrentDirectory, Path.GetFileNameWithoutExtension(genomeFasta) + "*"))
            {
                File.Delete(file);
            }
            foreach (var file in Directory.GetFiles(Environment.CurrentDirectory, Path.GetFileNameWithoutExtension(geneModelFile) + "*"))
            {
                File.Delete(file);
            }
            foreach (var file in Directory.GetFiles(Environment.CurrentDirectory, Path.GetFileNameWithoutExtension(proteinFasta) + "*"))
            {
                File.Delete(file);
            }
        }
    }
}