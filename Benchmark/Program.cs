﻿using Proteogenomics;
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
            SameProteins();
        }

        private static void SameProteins()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            // download and decompress references
            using (WebClient Client = new WebClient())
            {
                Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz");
                Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz", "Homo_sapiens.GRCh38.81.gff3.gz");
                Client.DownloadFile(@"ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz", "Homo_sapiens.GRCh38.pep.all.fa.gz");
            }
            string genomeFasta = "Homo_sapiens.GRCh38.dna.primary_assembly.fa";
            string geneModelFile = "Homo_sapiens.GRCh38.81.gff3";
            string proteinFasta = "Homo_sapiens.GRCh38.pep.all.fa";
            if (!File.Exists(genomeFasta))
            {
                using (FileStream stream = new FileStream("Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", FileMode.Open))
                using (GZipStream gunzip = new GZipStream(stream, CompressionMode.Decompress))
                using (var f = File.Create(proteinFasta))
                {
                    gunzip.CopyTo(f);
                }
            }
            if (!File.Exists(geneModelFile))
            {
                using (FileStream stream = new FileStream("Homo_sapiens.GRCh38.81.gff3.gz", FileMode.Open))
                using (GZipStream gunzip = new GZipStream(stream, CompressionMode.Decompress))
                using (var f = File.Create(proteinFasta))
                {
                    gunzip.CopyTo(f);
                }
            }
            if (!File.Exists(proteinFasta))
            {
                using (FileStream stream = new FileStream("Homo_sapiens.GRCh38.pep.all.fa.gz", FileMode.Open))
                using (GZipStream gunzip = new GZipStream(stream, CompressionMode.Decompress))
                using (var f = File.Create(proteinFasta))
                {
                    gunzip.CopyTo(f);
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