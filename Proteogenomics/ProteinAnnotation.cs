using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using System.Text.RegularExpressions;

namespace Proteogenomics
{
    public static class ProteinAnnotation
    {
        /// <summary>
        /// Transfers likely modifications from a list of proteins to another based on sequence similarity. Returns a list of new objects.
        /// </summary>
        /// <param name="destination"></param>
        /// <param name="source"></param>
        /// <returns></returns>
        public static List<Protein> TransferModifications(List<Protein> source, List<Protein> destination)
        {
            List<Protein> newProteins = new List<Protein>();
            Dictionary<string, List<Protein>> dictSource = ProteinDictionary(source);
            Dictionary<string, List<Protein>> dictDestination = ProteinDictionary(destination);

            List<string> seqsInCommon = dictDestination.Keys.Intersect(dictSource.Keys).ToList();
            List<string> destinationOnlySeqs = dictDestination.Keys.Except(dictSource.Keys).ToList();
            List<string> sourceOnlySeqs = dictSource.Keys.Except(dictDestination.Keys).ToList();

            foreach (string seq in seqsInCommon)
            {
                dictSource.TryGetValue(seq, out var sourceProteins);
                dictDestination.TryGetValue(seq, out var destinationProteins);
                if (sourceProteins == null || destinationProteins == null) { continue; } // shouldn't happen
                foreach (var ds in sourceProteins.SelectMany(s => s.DisulfideBonds)) { ds.Description = ds.Description ?? ""; } // shouldn't need this, and fixing in mzLib
                foreach (var destinationProtein in destinationProteins)
                {
                    newProteins.Add(
                        new Protein(

                            // keep these
                            seq,
                            destinationProtein.Accession,
                            name: destinationProtein.Name,
                            full_name: destinationProtein.FullName,
                            isDecoy: destinationProtein.IsDecoy,
                            isContaminant: destinationProtein.IsContaminant,
                            sequenceVariations: destinationProtein.SequenceVariations.ToList(),

                            // combine these
                            gene_names: destinationProtein.GeneNames.Concat(sourceProteins.SelectMany(p => p.GeneNames)).ToList(),

                            // transfer these
                            oneBasedModifications: CollapseMods(sourceProteins),
                            proteolysisProducts: new HashSet<ProteolysisProduct>(sourceProteins.SelectMany(p => p.ProteolysisProducts)).ToList(),
                            databaseReferences: new HashSet<DatabaseReference>(sourceProteins.SelectMany(p => p.DatabaseReferences)).ToList(),
                            disulfideBonds: new HashSet<DisulfideBond>(sourceProteins.SelectMany(p => p.DisulfideBonds)).ToList()
                        ));
                }
            }

            return newProteins.Concat(destinationOnlySeqs.SelectMany(s => dictDestination[s])).ToList();
        }

        public static List<ModificationWithLocation> GetUniProtMods(string spritzDirectory)
        {
            Loaders.LoadElements(Path.Combine(spritzDirectory, "elements.dat"));
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(spritzDirectory, "PSI-MOD.obo.xml"));
            return Loaders.LoadUniprot(Path.Combine(spritzDirectory, "ptmlist.txt"), Loaders.GetFormalChargesDictionary(psiModDeserialized)).ToList();
        }

        /// <summary>
        /// Collapses modification dictionaries from multiple proteins into one
        /// </summary>
        /// <param name="proteins"></param>
        /// <returns></returns>
        private static Dictionary<int, List<Modification>> CollapseMods(List<Protein> proteins)
        {
            Dictionary<int, HashSet<Modification>> result = new Dictionary<int, HashSet<Modification>>();
            foreach (var kv in proteins.SelectMany(p => p.OneBasedPossibleLocalizedModifications).ToList())
            {
                if (result.TryGetValue(kv.Key, out var modifications))
                {
                    foreach (var mod in kv.Value) { modifications.Add(mod); }
                }
                else
                {
                    result[kv.Key] = new HashSet<Modification>(kv.Value);
                }
            }
            return result.ToDictionary(kv => kv.Key, kv => kv.Value.ToList());
        }

        /// <summary>
        /// Creates a sequence-to-proteins dictionary
        /// </summary>
        /// <param name="proteins"></param>
        /// <returns></returns>
        private static Dictionary<string, List<Protein>> ProteinDictionary(List<Protein> proteins)
        {
            Dictionary<string, List<Protein>> dict = new Dictionary<string, List<Protein>>();
            foreach (Protein p in proteins)
            {
                if (dict.TryGetValue(p.BaseSequence, out var prots)) { prots.Add(p); }
                else { dict[p.BaseSequence] = new List<Protein> { p }; }
            }
            return dict;
        }

        /// <summary>
        /// Ensembl coding domain sequences (CDS) sometimes don't have start or stop codons annotated.
        /// The only way I can figure out how to tell which they are is to read in the protein FASTA and find the ones starting with X's or containing a stop codon '*'
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="proteinFastaPath"></param>
        /// <returns></returns>
        public static void GetImportantProteinAccessions(string proteinFastaPath, out Dictionary<string, string> proteinAccessionSequence, out HashSet<string> badProteinAccessions, 
            out Dictionary<string, string> selenocysteineProteinAccessions)
        {
            Regex transcriptAccession = new Regex(@"(transcript:)([A-Za-z0-9_.]+)"); // need to include transcript accessions for when a GTF file is used and transcript IDs become the protein IDs
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(proteinFastaPath, true, DecoyType.None, false,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblGeneNameRegex, null, out List<string> errors);
            proteinAccessionSequence = proteins.Select(p => new KeyValuePair<string, string>(p.Accession, p.BaseSequence))
                .Concat(proteins.Select(p => new KeyValuePair<string, string>(transcriptAccession.Match(p.FullName).Groups[2].Value, p.BaseSequence)))
                .ToDictionary(kv => kv.Key, kv => kv.Value);
            HashSet<string> badOnes = new HashSet<string>(proteins.Where(p => p.BaseSequence.Contains('X') || p.BaseSequence.Contains('*'))
                .SelectMany(p => new string[] { p.Accession, transcriptAccession.Match(p.FullName).Groups[2].Value }));
            badProteinAccessions = badOnes;
            selenocysteineProteinAccessions = proteins.Where(p => !badOnes.Contains(p.Accession) && p.BaseSequence.Contains('U')).ToDictionary(p => p.Accession, p => p.BaseSequence);
        }
    }
}