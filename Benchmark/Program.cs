using System;
using Bio;
using Proteomics;
using Proteogenomics;
using System.Collections.Generic;
using UsefulProteomicsDatabases;
using System.Linq;
using System.Diagnostics;

namespace Benchmark
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
        }

        public void SameProteins()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            // download and decompress references
            EnsemblDownloadsWrapper downloads = new EnsemblDownloadsWrapper();
            downloads.DownloadReferences(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch38");
            ProteinAnnotation.GetImportantProteinAccessions(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh38ProteinFastaFilename));

            Genome genome = new Genome(downloads.GenomeFastaPath);
            GeneModel geneModel = new GeneModel(genome, downloads.Gff3GeneModelPath);
            List<Protein> geneBasedProteins = geneModel.Translate(true, downloads.BadProteinAccessions, downloads.SelenocysteineProteinAccessions);
            List<Protein> pepDotAll = ProteinDbLoader.LoadProteinFasta(downloads.ProteinFastaPath, true, DecoyType.None, false,
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

            Console.WriteLine("Finished checking that all proteins are the same.");
            Console.WriteLine("Time elapsed: " + stopwatch.Elapsed.TotalHours.ToString() + ":" + stopwatch.Elapsed.Minutes.ToString() + ":" + stopwatch.Elapsed.Seconds.ToString());
            Console.WriteLine("Result: all proteins are " + (allAreEqual ? "" : "not") + " equal ");
        }
    }
}
