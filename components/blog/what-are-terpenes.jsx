import Link from 'next/link';

export default function WhatAreTerpenesContent() {
  return (
    <div className="prose prose-lg max-w-none">
      <p className="text-xl text-gray-700 leading-relaxed">
        That sharp, clean scent after rain in a pine forest? The bright burst of citrus when you peel an orange? The calming aroma of lavender sachets in your drawer? All terpenes. These molecules are responsible for most of the smells you encounter in nature—from the mint in your toothpaste to the hops in your beer—yet most people have never heard of them. If you&apos;ve seen &quot;terpenes&quot; mentioned on a cannabis label, essential oil bottle, or ingredient list and wondered what they actually <em>are</em>, here&apos;s the chemistry.
      </p>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Quick Answer: What Are Terpenes?</h2>

      <div className="bg-gray-50 border-l-4 border-gray-900 p-6 my-8">
        <p className="text-gray-700">
          Terpenes are volatile organic compounds built from repeating isoprene units (C₅H₈). They&apos;re the primary constituents of essential oils and responsible for the characteristic smells of most plants. Structurally, they&apos;re hydrocarbons—made only of carbon and hydrogen—classified by how many isoprene units they contain: monoterpenes (two units), sesquiterpenes (three units), and so on.
        </p>
      </div>

      <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-4">Key Takeaways</h3>
      <ul className="list-disc pl-6 mb-8 text-gray-700 space-y-2">
        <li>Terpenes are built from isoprene (C₅H₈) building blocks—not benzene rings, despite common confusion</li>
        <li>Over 30,000 terpene compounds exist in nature, making them the largest class of natural products</li>
        <li>The 3D shape of a terpene determines its smell—mirror-image molecules can smell completely different</li>
        <li>Most &quot;terpenes&quot; in commercial products are technically terpenoids (chemically modified terpenes)</li>
        <li>GC-MS (gas chromatography-mass spectrometry) is the gold standard for identifying terpenes in mixtures</li>
      </ul>

      <figure className="my-8">
        <img
          src="/terpene-isoprene.svg"
          alt="Isoprene molecular structure showing the basic C₅H₈ building block of all terpenes"
          title="Isoprene (C₅H₈) - the fundamental building block of terpenes"
          className="w-full max-w-sm mx-auto"
          loading="lazy"
          width="400"
          height="300"
        />
        <figcaption className="text-center text-sm text-gray-600 mt-2">
          Isoprene (C₅H₈) - the basic building block of all terpenes
        </figcaption>
      </figure>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">The Chemistry: What Terpenes Actually Are</h2>

      <p className="text-gray-700 leading-relaxed mb-4">
        Here&apos;s where most explanations get it wrong: terpenes are <strong>not</strong> aromatic compounds. They&apos;re not built around benzene rings, despite what the word &quot;aromatic&quot; might suggest. Instead, they&apos;re constructed from a simple five-carbon building block called isoprene.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        Isoprene has the molecular formula C₅H₈—five carbons arranged in a branched structure with a double bond. Think of it as a Lego brick. Connect two isoprene units together, and you get a monoterpene (C₁₀H₁₆). Three units? That&apos;s a sesquiterpene (C₁₅H₂₄). Four gives you a diterpene (C₂₀H₃₂). The pattern continues: triterpenes, tetraterpenes, polyterpenes.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        This modular structure is why terpenes are so diverse. The same building blocks can connect in different configurations, creating wildly different molecules. It&apos;s the chemical equivalent of how the same 26 letters create different words—except instead of meaning, you get different smells.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        The key property that makes terpenes smell: they&apos;re <strong>volatile</strong>. Their relatively small size and molecular structure mean they evaporate easily at room temperature. When you crush a basil leaf or peel a lemon, you&apos;re rupturing cells and releasing terpenes into the air, where they can travel to your nose.
      </p>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Common Terpenes and Where You&apos;ll Find Them</h2>

      <p className="text-gray-700 leading-relaxed mb-4">
        Over 30,000 terpenes have been identified in nature, but a handful dominate the smells you encounter daily. Here&apos;s what you&apos;re actually smelling:
      </p>

      <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Limonene</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Found in:</strong> Citrus peels (oranges, lemons, limes), caraway, dill<br/>
        <strong>Smell:</strong> Fresh, citrusy, slightly sweet<br/>
        <strong>Structure:</strong> Monoterpene (C₁₀H₁₆) with a six-membered ring
      </p>
      <p className="text-gray-700 leading-relaxed mb-4">
        Limonene makes up 90-95% of orange peel oil. It&apos;s the reason citrus-scented cleaning products actually smell like citrus. The compound is also used as a biodegradable solvent—it dissolves sticky residues and oils effectively, which is why it appears in industrial cleaners and adhesive removers.
      </p>

      <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Pinene</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Found in:</strong> Pine needles, rosemary, basil, cannabis<br/>
        <strong>Smell:</strong> Sharp, fresh, pine-like<br/>
        <strong>Structure:</strong> Monoterpene with a four-membered ring (strained, reactive)
      </p>
      <p className="text-gray-700 leading-relaxed mb-4">
        Pinene exists as two forms: α-pinene and β-pinene, which differ in where a double bond is located. Both smell similar—that characteristic pine scent—but α-pinene is more abundant. It&apos;s the most common terpene in nature globally, making up a significant portion of turpentine (yes, the paint thinner—it&apos;s distilled from pine resin).
      </p>

      <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Linalool</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Found in:</strong> Lavender, coriander, sweet basil, cinnamon<br/>
        <strong>Smell:</strong> Floral, slightly spicy, calming<br/>
        <strong>Structure:</strong> Monoterpene alcohol (technically a terpenoid)
      </p>
      <p className="text-gray-700 leading-relaxed mb-4">
        Linalool is why lavender is associated with relaxation—<a href="https://pubmed.ncbi.nlm.nih.gov/19962288/" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">it&apos;s been shown in studies to have mild sedative effects in mice</a>. It&apos;s also used extensively in perfumes and as a fragrance in soaps and detergents. Interestingly, <a href="https://pubmed.ncbi.nlm.nih.gov/26785373/" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">linalool oxidises when exposed to air, forming compounds that can cause skin sensitization</a>—which is why old lavender oil can be more irritating than fresh.
      </p>

      <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Myrcene</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Found in:</strong> Hops, mangoes, lemongrass, thyme<br/>
        <strong>Smell:</strong> Earthy, musky, slightly fruity<br/>
        <strong>Structure:</strong> Linear monoterpene (no ring structure)
      </p>
      <p className="text-gray-700 leading-relaxed mb-4">
        Myrcene is a major component of hop oil, contributing to beer&apos;s aroma. It&apos;s also why some cannabis strains smell earthy or musky. There&apos;s a persistent myth that eating mangoes before consuming cannabis enhances the effect because of myrcene content—this lacks scientific support but makes for good internet lore.
      </p>

      <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">β-Caryophyllene</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Found in:</strong> Black pepper, cloves, cannabis, hops<br/>
        <strong>Smell:</strong> Spicy, woody, peppery<br/>
        <strong>Structure:</strong> Sesquiterpene (C₁₅H₂₄) with a unique bicyclic structure
      </p>
      <p className="text-gray-700 leading-relaxed mb-4">
        This is what gives black pepper its spicy aroma. Unlike other terpenes, <a href="https://pubmed.ncbi.nlm.nih.gov/18574142/" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">β-caryophyllene can bind to cannabinoid receptors (CB2) in the human body</a>, making it technically a dietary cannabinoid—though it doesn&apos;t cause psychoactive effects. It&apos;s used in food flavouring and has been studied for potential anti-inflammatory properties.
      </p>

      <figure className="my-8">
        <img
          src="/terpene-structures-improved.png"
          alt="Molecular structures of isoprene, limonene, α-pinene, linalool, myrcene, and β-caryophyllene with stereochemistry shown"
          title="Common terpene molecular structures with stereochemistry: isoprene, limonene, α-pinene, linalool, myrcene, β-caryophyllene"
          className="w-full max-w-4xl mx-auto"
          loading="lazy"
          width="1500"
          height="800"
        />
        <figcaption className="text-center text-sm text-gray-600 mt-2">
          Common terpene structures showing different carbon skeletons built from isoprene units
        </figcaption>
      </figure>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Why Terpenes Smell: Molecular Shape and Receptors</h2>

      <p className="text-gray-700 leading-relaxed mb-4">
        Smell happens when volatile molecules bind to olfactory receptors in your nose—specialised proteins that change shape when the right molecule fits into them. Think of it like a lock and key: the molecule&apos;s 3D shape determines which receptors it activates.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        Here&apos;s where it gets fascinating: molecular shape is <em>everything</em>. Take limonene. It exists in two mirror-image forms called enantiomers—same atoms, same connections, but arranged as left-handed and right-handed versions. The right-handed form (d-limonene) smells like oranges. The left-handed form (l-limonene) smells like lemons or pine. Same molecular formula, completely different smell.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        This happens because your olfactory receptors are also chiral—they have handedness. A right-handed receptor can distinguish between right-handed and left-handed molecules, just like your right hand fits into a right-handed glove but not a left-handed one.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        <a href="https://www.nature.com/articles/nature12162" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">Humans have around 400 different olfactory receptor types</a>. Each terpene activates a different combination of these receptors, and your brain interprets that pattern as a specific smell. Change the molecular shape slightly—add a hydroxyl group, move a double bond—and you activate a different receptor combination, creating a different smell.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        This is why <a href="/blog/decaf-coffee-science" className="text-primary-600 hover:text-primary-700 underline">compounds extracted from coffee beans</a> can have such complex, varied aromas—many coffee volatiles are terpenoids, each contributing a slightly different note to the overall profile.
      </p>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Terpenes vs. Terpenoids: What&apos;s the Difference?</h2>

      <p className="text-gray-700 leading-relaxed mb-4">
        Technically, <strong>terpenes</strong> are pure hydrocarbons—only carbon and hydrogen. <strong>Terpenoids</strong> are terpenes that have been chemically modified, usually through oxidation, to include functional groups like alcohols (-OH), aldehydes (-CHO), ketones (=O), or esters (-COOR).
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        For example:
      </p>
      <ul className="list-disc pl-6 mb-4 text-gray-700">
        <li><strong>Limonene</strong> (pure hydrocarbon) → oxidises to <strong>carvone</strong> (contains a ketone group), which smells like spearmint or caraway</li>
        <li><strong>Menthene</strong> (pure hydrocarbon) → adds a hydroxyl group to become <strong>menthol</strong> (terpenoid), which has that cooling sensation</li>
        <li><strong>Geraniol</strong> (alcohol group present) is technically a terpenoid, even though it&apos;s often called a terpene</li>
      </ul>

      <p className="text-gray-700 leading-relaxed mb-4">
        In practice, most people use &quot;terpenes&quot; as a catch-all term. When a cannabis product lists &quot;terpene content,&quot; it&apos;s usually referring to both terpenes and terpenoids. The distinction matters more to chemists than to consumers—but if you want to be pedantic at parties, now you know.
      </p>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">How Scientists Identify Terpenes: GC-MS Analysis</h2>

      <p className="text-gray-700 leading-relaxed mb-4">
        When you see a cannabis product with a detailed terpene profile, or an essential oil labelled &quot;100% pure,&quot; that&apos;s been verified using <strong>gas chromatography-mass spectrometry (GC-MS)</strong>. It&apos;s the gold standard for analysing complex mixtures of volatile compounds.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        Here&apos;s how it works:
      </p>

      <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Step 1: Gas Chromatography (GC)</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        The sample—say, lavender oil—is vaporised and injected into a long, thin column (often 30+ meters coiled inside the instrument). The column is coated with a liquid stationary phase. As the vaporised sample travels through the column (carried by helium or nitrogen gas), different molecules move at different speeds.
      </p>
      <p className="text-gray-700 leading-relaxed mb-4">
        Small, volatile molecules like pinene travel quickly. Larger molecules like sesquiterpenes move slowly. Each compound exits the column at a different time—called its &quot;retention time.&quot; This separates the mixture into individual components.
      </p>

      <h3 className="text-2xl font-bold text-gray-900 mt-8 mb-3">Step 2: Mass Spectrometry (MS)</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        As each separated compound exits the column, it enters the mass spectrometer. Here, it&apos;s bombarded with electrons, fragmenting it into charged pieces. These fragments are separated by mass and detected, creating a unique pattern—a &quot;mass spectrum&quot;—like a molecular fingerprint.
      </p>
      <p className="text-gray-700 leading-relaxed mb-4">
        By comparing this fingerprint to a library of known spectra, the instrument identifies what the compound is. A typical GC-MS analysis of essential oil might separate 50-100+ individual compounds and identify most of them in a single run.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        This is how perfume chemists reverse-engineer competitor fragrances, how food scientists ensure flavour consistency, and how regulators verify that &quot;pure lavender oil&quot; isn&apos;t cut with synthetic linalool. It&apos;s also how cannabis testing labs quantify terpene content—though standardization in that industry remains inconsistent.
      </p>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Terpenes in Everyday Life</h2>

      <p className="text-gray-700 leading-relaxed mb-4">
        Once you know what terpenes are, you start seeing them everywhere:
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Perfume industry:</strong> Natural terpenes from flowers, woods, and resins form the backbone of many fragrances. Sandalwood oil (rich in sesquiterpenes) costs thousands per kilogram because it takes decades to grow sandalwood trees. Synthetic alternatives exist but don&apos;t quite match the complexity.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Cannabis:</strong> The &quot;entourage effect&quot; hypothesis suggests terpenes modulate how cannabinoids like THC and CBD affect the body. Scientific evidence remains limited, but terpene profiles definitely affect aroma and potentially user experience. Different strains have wildly different terpene compositions—some are pinene-heavy (alert, focused), others myrcene-dominant (sedating, relaxing).
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Essential oils:</strong> These are concentrated terpene mixtures extracted from plants. Tea tree oil, eucalyptus oil, peppermint oil—all predominantly terpenes and terpenoids. Their antimicrobial properties make them useful in natural cleaning products, though efficacy varies widely.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Food and beverage:</strong> That &quot;lemon fresh&quot; smell in dish soap? Limonene extracted from citrus peels. Hop aroma in beer? Myrcene, caryophyllene, and humulene. The cooling sensation in mint gum? Menthol, a terpenoid. The spice in ginger? Zingiberene and other sesquiterpenes. Many of these are also present in <a href="/blog/decaf-coffee-science" className="text-primary-600 hover:text-primary-700 underline">coffee&apos;s complex flavour profile</a>, contributing to the aromatic compounds that survive decaffeination.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        <strong>Cleaning products:</strong> Limonene is a powerful degreaser and solvent. It&apos;s biodegradable and less toxic than petroleum-based solvents, making it popular in &quot;green&quot; cleaners. It&apos;s also why orange peels can help remove sticky labels—rub citrus peel oil on adhesive residue, and it dissolves.
      </p>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Are Terpenes Safe?</h2>

      <p className="text-gray-700 leading-relaxed mb-4">
        Generally, yes—at the levels you encounter in food and natural products. <a href="https://www.fda.gov/food/food-ingredients-packaging/generally-recognized-safe-gras" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">Terpenes are &quot;Generally Recognized As Safe&quot; (GRAS) by the FDA</a> for use as food additives. You consume them every time you eat herbs, spices, or citrus.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        However, <strong>concentrated forms can be irritating</strong>. Essential oils are 50-100x more concentrated than the plant material they come from. <a href="https://www.poison.org/articles/essential-oils" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">Undiluted essential oils can cause skin irritation, allergic reactions, or respiratory issues in sensitive individuals</a>. <a href="https://www.aspca.org/pet-care/animal-poison-control/toxic-and-non-toxic-plants/tea-tree-oil" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">Pets, especially cats, metabolise terpenes poorly—tea tree oil and eucalyptus oil are particularly toxic to them</a>.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        There&apos;s also the &quot;natural = safe&quot; fallacy. Poison ivy produces urushiol (technically not a terpene, but a related compound). <a href="https://pubmed.ncbi.nlm.nih.gov/17448334/" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">Pennyroyal oil contains pulegone, a monoterpene that&apos;s hepatotoxic in high doses</a>. Thujone from wormwood is a neurotoxin. &quot;Natural&quot; doesn&apos;t mean harmless.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        For typical use—aromatherapy, cooking, cleaning products—terpenes are safe with proper dilution and ventilation. If you&apos;re using concentrated essential oils, follow dilution guidelines and keep them away from children and pets.
      </p>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">The Bottom Line</h2>

      <p className="text-gray-700 leading-relaxed mb-4">
        Terpenes are nature&apos;s vocabulary for smell. Built from simple five-carbon isoprene units, they create the vast majority of plant aromas you encounter—from the forest after rain to the basil in your pasta sauce. Their diversity comes from how those isoprene units connect, and their smell comes from how their 3D shape fits into your olfactory receptors.
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        Understanding terpenes gives you a new way to experience the world. That pine scent? You&apos;re smelling pinene evaporating from resin ducts in needles. That citrus burst? Limonene escaping from ruptured oil glands in the peel. That calming lavender? Linalool binding to receptors in your nose, triggering a cascade of neural signals your brain interprets as &quot;floral and relaxing.&quot;
      </p>

      <p className="text-gray-700 leading-relaxed mb-4">
        Chemistry doesn&apos;t diminish the experience—it deepens it.
      </p>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Frequently Asked Questions</h2>

      <div className="space-y-6 mb-12">
        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">What&apos;s the difference between terpenes and terpenoids?</h3>
          <p className="text-gray-700">
            Terpenes are hydrocarbons built from isoprene units—only carbon and hydrogen. Terpenoids are terpenes that have been chemically modified, usually through oxidation, adding functional groups like alcohols, aldehydes, or esters. For example, menthol is a terpenoid (it&apos;s modified from the terpene menthene by adding a hydroxyl group). In casual conversation, &quot;terpenes&quot; often refers to both.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">Are terpenes safe to inhale?</h3>
          <p className="text-gray-700">
            Most terpenes found naturally in foods and plants are safe at typical exposure levels. However, concentrated forms—like undiluted essential oils—can irritate airways or trigger allergic reactions in sensitive individuals. Proper ventilation is important when using products high in terpenes. Some people develop sensitivities to specific terpenes, particularly oxidised forms found in old essential oils.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">Why do terpenes smell?</h3>
          <p className="text-gray-700">
            Terpenes are volatile molecules that evaporate easily at room temperature and reach your nose. Once there, their specific three-dimensional shape determines how they bind to olfactory receptors—specialised proteins in your nasal cavity. Different molecular shapes activate different combinations of receptors, which your brain interprets as different smells. Even mirror-image versions of the same molecule (enantiomers) can smell completely different.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">What foods contain terpenes?</h3>
          <p className="text-gray-700">
            Nearly all aromatic plants contain terpenes. Common examples: citrus fruits (limonene), herbs like basil and rosemary (pinene, linalool), black pepper (caryophyllene), hops (myrcene), mangoes (myrcene), ginger (zingiberene), mint (menthol), cinnamon (cinnamaldehyde), and carrots (carotenes). Terpenes are responsible for most plant-based flavours and aromas in your kitchen.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">Do terpenes get you high?</h3>
          <p className="text-gray-700">
            No, terpenes themselves are not psychoactive. In cannabis, THC (a cannabinoid, not a terpene) causes the high. However, some research suggests terpenes may influence how cannabinoids affect the body—the &quot;entourage effect&quot;—though scientific evidence for this remains limited and debated. One exception: β-caryophyllene binds to cannabinoid receptors but doesn&apos;t produce psychoactive effects.
          </p>
        </div>
      </div>

      {/* Author byline */}
      <div className="mt-12 pt-8 border-t border-gray-200">
        <p className="text-gray-600 italic">
          Written by <strong>Miles Christou</strong>, chemistry graduate with expertise in analytical chemistry
        </p>
      </div>
    </div>
  );
}
