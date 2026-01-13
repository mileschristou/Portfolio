import { Link } from 'react-router-dom';
import SEO from '../../components/SEO';

const DecafCoffeeScience = () => {
  return (
    <>
      <SEO
        title="How is Coffee Decaffeinated? The Science Behind Your Decaf"
        description="Learn how coffee is decaffeinated using Swiss Water, CO₂, and solvent methods. Chemistry-backed guide to decaf processes, safety, and flavour."
      />

      {/* FAQ Schema for SEO */}
      <script type="application/ld+json">
        {JSON.stringify({
          "@context": "https://schema.org",
          "@type": "FAQPage",
          "mainEntity": [
            {
              "@type": "Question",
              "name": "Is decaf coffee completely caffeine-free?",
              "acceptedAnswer": {
                "@type": "Answer",
                "text": "No. Decaf coffee contains 2-7mg of caffeine per cup, compared to 70-140mg in regular coffee. FDA regulations allow up to 3% of original caffeine to remain."
              }
            },
            {
              "@type": "Question",
              "name": "Does decaffeination remove antioxidants from coffee?",
              "acceptedAnswer": {
                "@type": "Answer",
                "text": "Decaf retains most antioxidants, including chlorogenic acids and polyphenols. Studies show decaf still provides 70-80% of the antioxidant content of regular coffee."
              }
            },
            {
              "@type": "Question",
              "name": "Why does decaf coffee sometimes upset my stomach more than regular coffee?",
              "acceptedAnswer": {
                "@type": "Answer",
                "text": "The decaffeination process can slightly increase coffee's acidity and alter the balance of oils and compounds that affect digestion. Switching to Swiss Water or CO₂-processed decaf often helps."
              }
            },
            {
              "@type": "Question",
              "name": "Can pregnant women drink decaf coffee?",
              "acceptedAnswer": {
                "@type": "Answer",
                "text": "Most medical guidelines consider moderate decaf consumption (2-3 cups daily) safe during pregnancy. However, pregnant women should consult their healthcare provider about individual circumstances."
              }
            },
            {
              "@type": "Question",
              "name": "How can I tell which decaffeination method was used?",
              "acceptedAnswer": {
                "@type": "Answer",
                "text": "Packaging should specify the method, especially for Swiss Water or CO₂ processes. If unlabelled, it's likely a solvent-based method—usually methylene chloride in North America."
              }
            }
          ]
        })}
      </script>

      {/* Image Schema for SEO/AEO */}
      <script type="application/ld+json">
        {JSON.stringify([
          {
            "@context": "https://schema.org",
            "@type": "ImageObject",
            "contentUrl": "/caffeine-molecule.svg",
            "name": "Caffeine molecular structure",
            "description": "Caffeine (1,3,7-trimethylxanthine) molecular structure showing the xanthine core with three methyl groups",
            "encodingFormat": "image/svg+xml",
            "keywords": ["caffeine", "1,3,7-trimethylxanthine", "molecular structure", "decaffeination", "chemistry"]
          },
          {
            "@context": "https://schema.org",
            "@type": "ImageObject",
            "contentUrl": "/decaf-process-flow.svg",
            "name": "Decaffeination process flow diagram",
            "description": "Decaffeination process flow diagram showing green coffee beans, decaffeination, roasting, and final decaf coffee",
            "encodingFormat": "image/svg+xml",
            "keywords": ["decaffeination", "coffee", "process diagram", "decaf coffee"]
          },
          {
            "@context": "https://schema.org",
            "@type": "ImageObject",
            "contentUrl": "/swiss-water-process.svg",
            "name": "Swiss Water Process diagram",
            "description": "Swiss Water Process diagram showing caffeine molecules moving from coffee beans into caffeine-free water through osmosis",
            "encodingFormat": "image/svg+xml",
            "keywords": ["Swiss Water Process", "decaffeination", "osmosis", "caffeine extraction"]
          },
          {
            "@context": "https://schema.org",
            "@type": "ImageObject",
            "contentUrl": "/supercritical-co2.svg",
            "name": "Supercritical CO₂ extraction diagram",
            "description": "Supercritical CO₂ extraction process showing the continuous loop of CO₂ extracting caffeine from beans and being recycled",
            "encodingFormat": "image/svg+xml",
            "keywords": ["supercritical CO2", "decaffeination", "caffeine extraction", "CO2 process"]
          }
        ])}
      </script>

      <article className="py-12 px-4">
        <div className="max-w-3xl mx-auto">
          {/* Back link */}
          <Link to="/blog" className="inline-block text-gray-600 hover:text-gray-900 mb-8">
            ← Back to Blog
          </Link>

          {/* Article header */}
          <header className="mb-12">
            <h1 className="text-4xl md:text-5xl font-bold text-gray-900 mb-4 font-mono leading-tight">
              How is Coffee Decaffeinated? The Science Behind Your Decaf
            </h1>
            <div className="text-gray-600 text-sm">
              <time dateTime="2026-01-11">January 11, 2026</time>
              <span className="mx-2">•</span>
              <span>10 min read</span>
            </div>
          </header>

          {/* Article content */}
          <div className="prose prose-lg max-w-none">
            <p className="text-xl text-gray-700 leading-relaxed">
              Here's something that might surprise you: your "decaf" coffee isn't actually caffeine-free. By FDA standards, <a href="https://www.ncausa.org/Decaffeinated-Coffee" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">decaffeinated coffee can contain up to 3% of its original caffeine</a>—enough that a cup of decaf typically has 2-7mg of caffeine, compared to 70-140mg in regular coffee. But how do we remove 97% of the caffeine while keeping everything else that makes coffee taste like coffee? The answer involves some genuinely clever chemistry.
            </p>

            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Quick Answer: The Decaf Coffee Process Explained</h2>

            <div className="bg-gray-50 border-l-4 border-gray-900 p-6 my-8">
              <p className="font-semibold mb-2">How is coffee decaffeinated?</p>
              <p className="text-gray-700">
                Green (unroasted) coffee beans are treated with a solvent that selectively removes caffeine molecules. The main methods use organic solvents like ethyl acetate, supercritical carbon dioxide, or water-based processes that exploit caffeine's solubility. All methods remove caffeine before roasting, preserving most flavour compounds while extracting 96-99% of caffeine.
              </p>
            </div>

            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-4">Key Takeaways</h3>
            <ul className="list-disc pl-6 mb-8 text-gray-700 space-y-2">
              <li>Decaf coffee still contains 2-7mg of caffeine per cup (vs 70-140mg in regular coffee)</li>
              <li>Three main methods exist: solvent-based (cheapest), Swiss Water Process (chemical-free), and supercritical CO₂ (best flavour preservation)</li>
              <li>All FDA-approved decaffeination methods are safe—residual solvents evaporate during processing and roasting</li>
              <li>Decaffeination happens before roasting and inevitably removes some flavour compounds, which is why decaf tastes different</li>
              <li>Swiss Water and CO₂ methods cost 30-50% more but preserve significantly more flavour complexity</li>
            </ul>

            <figure className="my-8">
              <img
                src="/decaf-process-flow.svg"
                alt="Decaffeination process flow diagram showing green coffee beans, decaffeination, roasting, and final decaf coffee"
                title="Decaffeination process flow diagram"
                className="w-full max-w-3xl mx-auto"
                loading="lazy"
                width="1100"
                height="400"
              />
            </figure>

            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">A Brief (and Slightly Toxic) History</h2>

            <p className="text-gray-700 leading-relaxed mb-4">
              The first commercially successful decaffeination process was patented by German coffee merchant Ludwig Roselius in 1905. His breakthrough was realizing that caffeine could be extracted from steamed coffee beans. The catch? His solvent of choice was benzene—yes, the same benzene we now know as a carcinogen.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              Roselius stumbled onto his method after a shipment of coffee beans was soaked in seawater during transport, inadvertently removing much of the caffeine. He refined this into a process using benzene to dissolve caffeine from steamed beans. The resulting brand, Sanka (short for "sans caffeine"), dominated the decaf market for decades.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              By the 1970s, we'd figured out that regularly drinking benzene-treated coffee wasn't ideal, leading to <a href="https://cen.acs.org/food/food-science/coffee-decaffeinated-safe-drink/102/i27" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">the modern methods we use today</a>. Each exploits the same fundamental principle: caffeine is more soluble in certain solvents than the compounds that give coffee its flavour.
            </p>

            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Method 1: Solvent-Based Decaffeination</h2>

            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">The Chemistry of Selective Extraction</h3>

            <p className="text-gray-700 leading-relaxed mb-4">
              Why can we remove caffeine without stripping away everything else? It comes down to molecular properties. Caffeine (1,3,7-trimethylxanthine, if we're being formal) is a polar alkaloid that dissolves readily in certain organic solvents. The flavour compounds in coffee—oils, sugars, acids—have different polarities and molecular weights, making them less soluble in these same solvents.
            </p>

            <figure className="my-8">
              <img
                src="/caffeine-molecule.svg"
                alt="Caffeine molecular structure (1,3,7-trimethylxanthine) showing the xanthine core with three methyl groups"
                title="Caffeine (1,3,7-trimethylxanthine) molecular structure"
                className="w-full max-w-md mx-auto"
                loading="lazy"
                width="500"
                height="400"
              />
              <figcaption className="text-center text-sm text-gray-600 mt-2">
                Caffeine (1,3,7-trimethylxanthine) - the molecule we're trying to remove
              </figcaption>
            </figure>

            <p className="text-gray-700 leading-relaxed mb-4">
              Think of it like panning for gold: you're exploiting a property (density for gold, solubility for caffeine) that lets you separate what you want from everything else.
            </p>

            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Direct Solvent Method</h3>

            <p className="text-gray-700 leading-relaxed mb-4">
              In the direct method, green coffee beans are steamed to open their pores, then repeatedly rinsed with a solvent—typically methylene chloride or ethyl acetate. The caffeine dissolves into the solvent, which is drained away. After 8-12 hours of this treatment, the beans are steamed again to evaporate any residual solvent.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>Methylene chloride</strong> is exceptionally good at dissolving caffeine while leaving flavour compounds alone. It's also the same chemical used in paint strippers, which understandably makes some people nervous. However, it has a boiling point of 40°C (104°F), meaning it evaporates completely during the subsequent steaming and roasting process. <a href="https://www.law.cornell.edu/cfr/text/21/173.255" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">FDA regulations require residual levels below 10 parts per million</a>—roasted decaf coffee typically contains less than 1ppm.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>Ethyl acetate</strong> offers a marketing advantage: it occurs naturally in fruits, allowing brands to advertise "naturally decaffeinated" coffee. Chemically, it works similarly to methylene chloride but is slightly less selective, potentially removing more flavour compounds along with the caffeine.
            </p>

            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Indirect Solvent Method</h3>

            <p className="text-gray-700 leading-relaxed mb-4">
              The indirect method is a clever workaround to avoid direct solvent contact with beans. Coffee beans are soaked in hot water, which extracts both caffeine and flavour compounds. The beans are removed, and the water is treated with solvent to extract caffeine. The now-caffeine-free but flavour-rich water is reintroduced to the beans, which reabsorb the flavour compounds.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              This method typically preserves more flavour than direct solvent extraction, but it's more time-consuming and expensive.
            </p>


            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Is Decaf Coffee Safe?</h3>

            <p className="text-gray-700 leading-relaxed mb-4">
              The residual solvent question comes up constantly. Here's the straightforward answer: modern solvent-based decaf is safe. Both methylene chloride and ethyl acetate are volatile (low boiling point), meaning they evaporate during processing. Roasting temperatures exceed 200°C (400°F)—well above what's needed to eliminate any traces. Multiple studies have found no detectable solvent residues in finished decaf coffee.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              The European Union, known for strict food safety standards, permits both methods. If there were legitimate health concerns, that wouldn't be the case.
            </p>

            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Method 2: The Swiss Water Process</h2>

            <p className="text-gray-700 leading-relaxed mb-4">
              The <a href="https://www.swisswater.com/pages/coffee-decaffeination-process" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">Swiss Water Process</a> is the darling of specialty coffee roasters, marketed as "chemical-free" decaffeination. That's technically true—no organic solvents involved—but it's still very much a chemical process. It just uses water and osmosis instead of industrial solvents.
            </p>

            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">How Water Removes Caffeine</h3>

            <p className="text-gray-700 leading-relaxed mb-4">
              Green coffee beans are soaked in hot water, which extracts caffeine along with flavour compounds (sugars, oils, acids). Here's the clever bit: that water is then passed through activated carbon filters with pores precisely sized to trap caffeine molecules while letting smaller flavour compounds pass through.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              Now you have caffeine-free water saturated with coffee flavour compounds—called "green coffee extract" or GCE. Fresh beans are soaked in this GCE. Because the water is already saturated with flavour compounds, osmosis doesn't pull them from the beans. But the water has no caffeine, creating a concentration gradient that pulls caffeine out of the beans and into the water.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              It's essentially controlled diffusion: caffeine molecules migrate from high concentration (inside beans) to low concentration (in the water) until equilibrium is reached. The process repeats with fresh GCE until 99% of caffeine is removed.
            </p>

            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Why Specialty Coffee Prefers It</h3>

            <p className="text-gray-700 leading-relaxed mb-4">
              The Swiss Water Process preserves flavour compounds better than solvent methods—particularly the delicate, complex notes that specialty coffee nerds care about. There's also the marketing advantage: "chemical-free" sounds better than "methylene chloride-treated," even though both are safe.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              The downside? Cost. It's significantly more expensive than solvent methods, which is why Swiss Water decaf often costs 30-50% more than conventional decaf.
            </p>

            <figure className="my-8">
              <img
                src="/swiss-water-process.svg"
                alt="Swiss Water Process diagram showing caffeine molecules moving from coffee beans into caffeine-free water through osmosis"
                title="Swiss Water Process - chemical-free decaffeination using osmosis"
                className="w-full max-w-3xl mx-auto"
                loading="lazy"
                width="1100"
                height="500"
              />
            </figure>

            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Method 3: Supercritical CO₂ Decaffeination</h2>

            <p className="text-gray-700 leading-relaxed mb-4">
              This is where chemistry gets genuinely fascinating. Carbon dioxide—the gas you exhale—becomes a powerful, selective solvent under the right conditions.
            </p>

            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">What "Supercritical" Means</h3>

            <p className="text-gray-700 leading-relaxed mb-4">
              Above a certain pressure and temperature (73 atmospheres and 31°C for CO₂), carbon dioxide enters a supercritical state: it's not quite a gas, not quite a liquid. It has the density of a liquid (so it can dissolve things) but the viscosity of a gas (so it penetrates materials easily).
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              Supercritical CO₂ is highly selective for caffeine. Beans are placed in a high-pressure vessel, supercritical CO₂ is circulated through them, and caffeine dissolves into the CO₂. The caffeine-laden CO₂ is then passed through water or activated carbon to separate the caffeine, and the CO₂ is recycled back through the beans.
            </p>

            <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Advantages and Limitations</h3>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>Pros:</strong>
            </p>
            <ul className="list-disc pl-6 mb-4 text-gray-700">
              <li>Extremely selective—removes caffeine while preserving flavour better than any other method</li>
              <li>CO₂ is non-toxic, naturally occurring, and fully recyclable</li>
              <li>Environmentally friendly—no chemical waste</li>
            </ul>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>Cons:</strong>
            </p>
            <ul className="list-disc pl-6 mb-4 text-gray-700">
              <li>Requires expensive high-pressure equipment (think industrial-scale pressure cookers)</li>
              <li>Higher operating costs than solvent methods</li>
              <li>Only economical at large scale</li>
            </ul>

            <p className="text-gray-700 leading-relaxed mb-4">
              This method is commonly used for high-quality decaf and is popular in Europe, where regulations favor it over solvent-based methods.
            </p>

            <figure className="my-8">
              <img
                src="/supercritical-co2.svg"
                alt="Supercritical CO₂ extraction process showing the continuous loop of CO₂ extracting caffeine from beans and being recycled"
                title="Supercritical CO₂ extraction - the most selective decaffeination method"
                className="w-full max-w-3xl mx-auto"
                loading="lazy"
                width="1100"
                height="550"
              />
            </figure>

            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Why Decaf Tastes Different</h2>

            <p className="text-gray-700 leading-relaxed mb-4">
              Even with perfect decaffeination, decaf never tastes quite like regular coffee. Why?
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>Flavor compound loss:</strong> No method is perfectly selective. Some volatile flavour compounds inevitably get removed along with caffeine, particularly lighter, fruity notes. The Swiss Water and supercritical CO₂ methods minimise this, but some loss is unavoidable.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>Roasting differences:</strong> Decaffeination makes beans more porous and brittle, changing how they respond to roasting. Roasters often need to adjust time and temperature, which affects flavour development. Beans that have been through decaffeination roast faster and less evenly.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>Missing caffeine's bitterness:</strong> Caffeine itself contributes to coffee's characteristic bitterness. Remove it, and the flavour profile shifts—often perceived as flatter or less complex.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>Bean quality:</strong> Historically, lower-quality beans were used for decaf (why waste premium beans on decaf?). That's changing with specialty decaf, but it's still a factor in commodity-grade decaf.
            </p>

            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Which Decaffeination Method is Best?</h2>

            <p className="text-gray-700 leading-relaxed mb-4">
              There's no single "best" method—it depends on priorities:
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>For flavour:</strong> Supercritical CO₂ or Swiss Water Process preserve the most complexity. Expect to pay premium prices.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>For cost:</strong> Solvent-based methods (especially direct methylene chloride) are cheapest and most widely available. Flavor loss is noticeable but acceptable for everyday drinking.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>For "natural" preferences:</strong> Swiss Water Process is the only method that doesn't use organic solvents, though "chemical-free" is misleading (everything is chemicals, including water).
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>For environmental impact:</strong> Supercritical CO₂ wins—it's recyclable and produces no chemical waste. Swiss Water is water-intensive. Solvent methods generate chemical waste that requires disposal.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              <strong>For safety:</strong> All methods are safe when properly executed. Residual solvents in finished coffee are negligible to nonexistent.
            </p>


            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">The Bottom Line</h2>

            <p className="text-gray-700 leading-relaxed mb-4">
              Coffee decaffeination is a surprisingly elegant application of chemistry: exploiting solubility differences to remove one specific compound while preserving hundreds of others. Whether through solvents, water, or supercritical fluids, the goal is the same—give you the ritual and flavour of coffee without the 2am wide-awake staring at the ceiling.
            </p>

            <p className="text-gray-700 leading-relaxed mb-4">
              Modern decaf is safer, better-tasting, and more environmentally friendly than ever. And while it won't perfectly replicate regular coffee, the gap is narrowing—especially if you're willing to pay for Swiss Water or CO₂-processed beans from a quality roaster.
            </p>

            <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Frequently Asked Questions</h2>

            <div className="space-y-6 mb-12">
              <div className="border-l-4 border-gray-400 pl-6">
                <h3 className="text-xl font-bold text-gray-900 mb-2">Is decaf coffee completely caffeine-free?</h3>
                <p className="text-gray-700">
                  No. Decaf coffee contains 2-7mg of caffeine per cup, compared to 70-140mg in regular coffee. FDA regulations allow up to 3% of original caffeine to remain. For most people, this trace amount won't cause noticeable effects, but those extremely sensitive to caffeine should be aware.
                </p>
              </div>

              <div className="border-l-4 border-gray-400 pl-6">
                <h3 className="text-xl font-bold text-gray-900 mb-2">Does decaffeination remove antioxidants from coffee?</h3>
                <p className="text-gray-700">
                  Decaf retains most antioxidants, including chlorogenic acids and polyphenols. Some loss occurs during processing, but <a href="https://www.ncausa.org/Decaffeinated-Coffee" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">studies show decaf still provides 70-80% of the antioxidant content</a> of regular coffee. The health benefits are largely preserved.
                </p>
              </div>

              <div className="border-l-4 border-gray-400 pl-6">
                <h3 className="text-xl font-bold text-gray-900 mb-2">Why does decaf coffee sometimes upset my stomach more than regular coffee?</h3>
                <p className="text-gray-700">
                  The decaffeination process can slightly increase coffee's acidity and alter the balance of oils and compounds that affect digestion. Additionally, some solvent-based methods may leave behind trace compounds that irritate sensitive stomachs. Switching to Swiss Water or CO₂-processed decaf often helps.
                </p>
              </div>

              <div className="border-l-4 border-gray-400 pl-6">
                <h3 className="text-xl font-bold text-gray-900 mb-2">Can pregnant women drink decaf coffee?</h3>
                <p className="text-gray-700">
                  Most medical guidelines consider moderate decaf consumption (2-3 cups daily) safe during pregnancy. The minimal caffeine content falls well below recommended limits. However, pregnant women should consult their healthcare provider about individual circumstances and overall caffeine intake from all sources.
                </p>
              </div>

              <div className="border-l-4 border-gray-400 pl-6">
                <h3 className="text-xl font-bold text-gray-900 mb-2">How can I tell which decaffeination method was used?</h3>
                <p className="text-gray-700">
                  Packaging should specify the method, especially for Swiss Water or CO₂ processes (premium methods that brands advertise). If unlabelled, it's likely a solvent-based method—usually methylene chloride in North America, ethyl acetate if marketed as "naturally decaffeinated." When in doubt, contact the roaster directly.
                </p>
              </div>
            </div>

            {/* Author byline */}
            <div className="mt-12 pt-8 border-t border-gray-200">
              <p className="text-gray-600 italic">
                Written by <strong>Miles Christou</strong>, chemistry graduate and coffee lover
              </p>
            </div>
          </div>

          {/* Back to blog link */}
          <div className="mt-8 pt-8 border-t border-gray-200">
            <Link to="/blog" className="inline-block text-gray-900 font-medium hover:underline">
              ← Back to all posts
            </Link>
          </div>
        </div>
      </article>
    </>
  );
};

export default DecafCoffeeScience;
