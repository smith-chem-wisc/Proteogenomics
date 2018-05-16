using NUnit.Framework;
using Proteogenomics;
using System.IO;

namespace TestProteogenomics
{
    [TestFixture]
    public class Test
    {
        [Test]
        public void CountReads()
        {
            Assert.AreEqual(3970, new FastqProperties(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper.fastq")).ReadCount);
        }
    }
}