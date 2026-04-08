export default function AnalyticalMethodsContent() {
  return (
    <div className="prose prose-lg max-w-none">
      <p className="text-xl text-gray-700 leading-relaxed">
        Walk into any modern analytical QC laboratory and you&apos;ll find three instruments that together can identify and quantify almost anything in a sample — metals at parts-per-trillion, ions at parts-per-million, and volatile organic compounds across a huge range of concentrations. ICP, IC, and GC are the workhorses of environmental and industrial QC. This guide explains how each one actually works, from first principles.
      </p>

      {/* ===================== QUICK SUMMARY ===================== */}
      <div className="bg-gray-50 border-l-4 border-gray-900 p-6 my-8">
        <p className="text-gray-800 font-semibold mb-2">Quick Reference</p>
        <ul className="text-gray-700 space-y-1 mb-0">
          <li><strong>ICP</strong> (Inductively Coupled Plasma) — elemental analysis; detects metals and trace elements in solution</li>
          <li><strong>IC</strong> (Ion Chromatography) — separates and quantifies dissolved ions; great for anions and cations in water</li>
          <li><strong>GC</strong> (Gas Chromatography) — separates volatile organic compounds; paired with FID or MS detection</li>
        </ul>
      </div>

      {/* ===================== ICP ===================== */}
      <h2 className="text-3xl font-bold text-gray-900 mt-14 mb-4 font-mono">ICP — Inductively Coupled Plasma</h2>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">What is it for?</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        ICP is used for <strong>elemental analysis</strong> — determining which metals and other elements are present in a sample, and at what concentration. It can detect dozens of elements simultaneously, from major constituents down to trace levels in the parts-per-billion (μg/L) or even parts-per-trillion (ng/L) range. Common applications include metals in drinking water, heavy metals in soil, elemental composition of industrial products, and nutrient analysis in fertilisers (phosphorus, potassium, calcium, magnesium, etc.).
      </p>
      <p className="text-gray-700 leading-relaxed mb-4">
        There are two main variants: <strong>ICP-OES</strong> (Optical Emission Spectrometry, also called ICP-AES) and <strong>ICP-MS</strong> (Mass Spectrometry). ICP-OES is more common in routine QC labs and handles most environmental and industrial work. ICP-MS goes to lower detection limits and is used for ultra-trace analysis.
      </p>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">The core idea: atomise everything, then identify the atoms</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        The fundamental challenge of elemental analysis is that you need to break every molecule in the sample apart into individual atoms, get those atoms excited, and then measure what they emit (OES) or weigh them (MS). ICP achieves the first two steps with an argon plasma — one of the hottest analytical tools in routine chemistry.
      </p>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">How it works, step by step</h3>

      {/* ICP Diagram */}
      <figure className="my-8">
        <svg viewBox="0 0 820 320" className="w-full max-w-3xl mx-auto" aria-label="ICP-OES instrument diagram">
          {/* Background */}
          <rect width="820" height="320" fill="#f9fafb" rx="8"/>

          {/* --- NEBULISER --- */}
          <rect x="30" y="120" width="100" height="80" rx="6" fill="#e5e7eb" stroke="#374151" strokeWidth="1.5"/>
          <text x="80" y="152" textAnchor="middle" fontSize="12" fontWeight="bold" fill="#111827">Nebuliser</text>
          <text x="80" y="168" textAnchor="middle" fontSize="10" fill="#374151">+ Spray</text>
          <text x="80" y="181" textAnchor="middle" fontSize="10" fill="#374151">Chamber</text>
          {/* sample input arrow */}
          <line x1="80" y1="120" x2="80" y2="95" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr)"/>
          <text x="80" y="88" textAnchor="middle" fontSize="10" fill="#374151">Liquid sample</text>
          {/* Ar gas */}
          <line x1="30" y1="160" x2="10" y2="160" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr)"/>
          <text x="8" y="155" textAnchor="end" fontSize="9" fill="#374151">Ar</text>

          {/* arrow to torch */}
          <line x1="130" y1="160" x2="175" y2="160" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr)"/>
          <text x="152" y="153" textAnchor="middle" fontSize="9" fill="#374151">aerosol</text>

          {/* --- TORCH + PLASMA --- */}
          {/* torch body */}
          <rect x="175" y="90" width="70" height="140" rx="6" fill="#dbeafe" stroke="#1d4ed8" strokeWidth="1.5"/>
          <text x="210" y="115" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#1e40af">Torch</text>
          {/* plasma flame */}
          <ellipse cx="210" cy="195" rx="22" ry="32" fill="#fef08a" stroke="#ca8a04" strokeWidth="1.5" opacity="0.9"/>
          <ellipse cx="210" cy="185" rx="14" ry="22" fill="#fde68a" stroke="#b45309" strokeWidth="1" opacity="0.8"/>
          <ellipse cx="210" cy="178" rx="7" ry="14" fill="#fed7aa" stroke="#ea580c" strokeWidth="1" opacity="0.9"/>
          <text x="210" y="240" textAnchor="middle" fontSize="10" fontWeight="bold" fill="#92400e">~8000 K</text>
          <text x="210" y="253" textAnchor="middle" fontSize="9" fill="#92400e">plasma</text>
          {/* RF coil indication */}
          <path d="M 175 165 Q 188 157 200 165 Q 213 173 225 165 Q 238 157 245 165" fill="none" stroke="#b45309" strokeWidth="2" strokeDasharray="3,2"/>
          <text x="285" y="168" textAnchor="start" fontSize="9" fill="#b45309">RF coil</text>
          <line x1="245" y1="165" x2="280" y2="165" stroke="#b45309" strokeWidth="1" strokeDasharray="2,2"/>

          {/* arrow to spectrometer */}
          <line x1="245" y1="160" x2="310" y2="160" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr)"/>
          <text x="277" y="150" textAnchor="middle" fontSize="9" fill="#374151">light / ions</text>

          {/* --- OES SPECTROMETER --- */}
          <rect x="310" y="100" width="160" height="120" rx="6" fill="#f3e8ff" stroke="#7e22ce" strokeWidth="1.5"/>
          <text x="390" y="128" textAnchor="middle" fontSize="12" fontWeight="bold" fill="#6b21a8">Spectrometer</text>
          <text x="390" y="145" textAnchor="middle" fontSize="10" fill="#6b21a8">(OES variant)</text>
          {/* grating lines */}
          <line x1="330" y1="165" x2="470" y2="165" stroke="#a855f7" strokeWidth="1"/>
          <line x1="340" y1="155" x2="340" y2="195" stroke="#a855f7" strokeWidth="1.5"/>
          <line x1="360" y1="155" x2="360" y2="195" stroke="#ef4444" strokeWidth="1.5"/>
          <line x1="380" y1="155" x2="380" y2="195" stroke="#22c55e" strokeWidth="1.5"/>
          <line x1="400" y1="155" x2="400" y2="195" stroke="#3b82f6" strokeWidth="1.5"/>
          <line x1="420" y1="155" x2="420" y2="195" stroke="#f59e0b" strokeWidth="1.5"/>
          <text x="390" y="212" textAnchor="middle" fontSize="9" fill="#7e22ce">diffraction grating → CCD</text>

          {/* arrow to detector */}
          <line x1="470" y1="160" x2="540" y2="160" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr)"/>

          {/* --- DETECTOR / DATA --- */}
          <rect x="540" y="115" width="120" height="90" rx="6" fill="#dcfce7" stroke="#166534" strokeWidth="1.5"/>
          <text x="600" y="147" textAnchor="middle" fontSize="12" fontWeight="bold" fill="#166534">Detector</text>
          <text x="600" y="163" textAnchor="middle" fontSize="10" fill="#166534">Wavelength →</text>
          <text x="600" y="178" textAnchor="middle" fontSize="10" fill="#166534">Element ID</text>
          <text x="600" y="193" textAnchor="middle" fontSize="10" fill="#166534">Intensity → conc.</text>

          {/* arrow to result */}
          <line x1="660" y1="160" x2="720" y2="160" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr)"/>
          <rect x="720" y="130" width="80" height="60" rx="6" fill="#fff7ed" stroke="#c2410c" strokeWidth="1.5"/>
          <text x="760" y="158" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#c2410c">Result</text>
          <text x="760" y="172" textAnchor="middle" fontSize="9" fill="#c2410c">e.g. Pb: 2.3 μg/L</text>

          {/* Labels underneath */}
          <text x="80" y="290" textAnchor="middle" fontSize="10" fill="#6b7280">① Convert to</text>
          <text x="80" y="302" textAnchor="middle" fontSize="10" fill="#6b7280">fine aerosol</text>
          <text x="210" y="290" textAnchor="middle" fontSize="10" fill="#6b7280">② Atomise &amp;</text>
          <text x="210" y="302" textAnchor="middle" fontSize="10" fill="#6b7280">excite atoms</text>
          <text x="390" y="290" textAnchor="middle" fontSize="10" fill="#6b7280">③ Separate by</text>
          <text x="390" y="302" textAnchor="middle" fontSize="10" fill="#6b7280">wavelength</text>
          <text x="600" y="290" textAnchor="middle" fontSize="10" fill="#6b7280">④ Identify &amp;</text>
          <text x="600" y="302" textAnchor="middle" fontSize="10" fill="#6b7280">quantify</text>

          {/* arrowhead marker */}
          <defs>
            <marker id="arr" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
              <path d="M0,0 L0,6 L8,3 z" fill="#374151"/>
            </marker>
          </defs>
        </svg>
        <figcaption className="text-center text-sm text-gray-600 mt-2">
          ICP-OES schematic: liquid sample → nebuliser → argon plasma torch → spectrometer → result
        </figcaption>
      </figure>

      <ol className="list-decimal pl-6 mb-6 text-gray-700 space-y-3">
        <li>
          <strong>Sample introduction.</strong> The liquid sample (usually an acidified aqueous solution) is pumped into a <em>nebuliser</em>, which converts it into a fine aerosol mist using a flow of argon gas. A spray chamber removes the largest droplets — only the finest particles continue, because only they will fully vaporise in the plasma.
        </li>
        <li>
          <strong>The plasma torch.</strong> The aerosol is swept into a quartz torch surrounded by a copper RF (radio-frequency) coil. The coil oscillates at ~27 MHz, creating a rapidly changing magnetic field. This field induces currents in the ionised argon gas, heating it to between 6,000 and 10,000 K — hotter than the surface of the sun. At these temperatures, molecules are completely destroyed, and every element is atomised and then ionised.
        </li>
        <li>
          <strong>Emission (OES) or mass separation (MS).</strong> In ICP-OES, excited atoms emit light as electrons fall back to lower energy levels. Because each element has a unique set of emission wavelengths (a spectral fingerprint), a diffraction grating separates the light and a CCD detector records intensities at hundreds of wavelengths simultaneously. In ICP-MS, the plasma ions are instead extracted through a series of cones into a vacuum, where they are separated by their mass-to-charge ratio — giving even lower detection limits and isotopic information.
        </li>
        <li>
          <strong>Quantification.</strong> Concentration is determined by comparing signal intensity against calibration standards of known concentration. The relationship is linear over several orders of magnitude, which is one of ICP&apos;s great strengths.
        </li>
      </ol>

      <div className="bg-gray-50 border border-gray-200 p-5 my-6">
        <p className="text-gray-700 font-semibold mb-1">Key things to remember about ICP</p>
        <ul className="text-gray-700 space-y-1 text-sm">
          <li>— Samples must be in liquid form (solid samples are digested in acid first)</li>
          <li>— Measures elements, not molecules — it cannot tell you what compound an element was in</li>
          <li>— ICP-OES: detection limits typically μg/L (ppb) range; ICP-MS: ng/L (ppt) range</li>
          <li>— Internal standards (e.g. indium, rhodium) are added to correct for instrument drift and matrix effects</li>
          <li>— Matrix matching of calibrants to samples is important — dissolved solids and acids affect the signal</li>
        </ul>
      </div>

      {/* ===================== IC ===================== */}
      <h2 className="text-3xl font-bold text-gray-900 mt-14 mb-4 font-mono">IC — Ion Chromatography</h2>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">What is it for?</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        IC separates and quantifies dissolved <strong>ions</strong> in aqueous samples. In environmental and industrial QC labs it is the go-to technique for anions such as fluoride (F⁻), chloride (Cl⁻), nitrate (NO₃⁻), sulphate (SO₄²⁻), and phosphate (PO₄³⁻), and for cations like sodium (Na⁺), potassium (K⁺), calcium (Ca²⁺), and magnesium (Mg²⁺). For a company working with fertilisers or minerals, IC is particularly valuable — phosphate, nitrate, sulphate, and potassium are all routine targets.
      </p>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">The core idea: separate ions by how strongly they cling to a charged column</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        IC is a form of liquid chromatography. The sample is carried through a column packed with ion-exchange resin — a material bearing fixed charged groups on its surface. Ions in the sample compete with the eluent (the mobile phase, typically a dilute carbonate or hydroxide solution) for binding sites on the resin. Ions with a stronger attraction to the resin spend more time bound to it, travel more slowly through the column, and elute later. This differences in retention time is what separates them.
      </p>

      {/* IC Diagram */}
      <figure className="my-8">
        <svg viewBox="0 0 860 300" className="w-full max-w-3xl mx-auto" aria-label="Ion chromatography instrument diagram">
          <rect width="860" height="300" fill="#f9fafb" rx="8"/>

          <defs>
            <marker id="arr2" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
              <path d="M0,0 L0,6 L8,3 z" fill="#374151"/>
            </marker>
          </defs>

          {/* Eluent reservoir */}
          <rect x="20" y="90" width="80" height="80" rx="6" fill="#dbeafe" stroke="#1d4ed8" strokeWidth="1.5"/>
          <text x="60" y="124" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#1e40af">Eluent</text>
          <text x="60" y="140" textAnchor="middle" fontSize="9" fill="#1e40af">reservoir</text>
          <text x="60" y="155" textAnchor="middle" fontSize="9" fill="#1e40af">(e.g. Na₂CO₃)</text>

          <line x1="100" y1="130" x2="140" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr2)"/>

          {/* Pump */}
          <rect x="140" y="105" width="70" height="50" rx="6" fill="#e5e7eb" stroke="#374151" strokeWidth="1.5"/>
          <text x="175" y="132" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#111827">Pump</text>
          <text x="175" y="146" textAnchor="middle" fontSize="9" fill="#374151">high pressure</text>

          <line x1="210" y1="130" x2="245" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr2)"/>

          {/* Injector */}
          <rect x="245" y="105" width="70" height="50" rx="6" fill="#fef9c3" stroke="#ca8a04" strokeWidth="1.5"/>
          <text x="280" y="128" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#854d0e">Sample</text>
          <text x="280" y="142" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#854d0e">Injector</text>
          {/* sample arrow in */}
          <line x1="280" y1="105" x2="280" y2="75" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr2)"/>
          <text x="280" y="68" textAnchor="middle" fontSize="9" fill="#374151">sample in</text>

          <line x1="315" y1="130" x2="355" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr2)"/>

          {/* Column */}
          <rect x="355" y="80" width="80" height="100" rx="8" fill="#f3e8ff" stroke="#7e22ce" strokeWidth="2"/>
          <text x="395" y="120" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#6b21a8">Column</text>
          <text x="395" y="136" textAnchor="middle" fontSize="9" fill="#6b21a8">ion-exchange</text>
          <text x="395" y="149" textAnchor="middle" fontSize="9" fill="#6b21a8">resin</text>
          {/* resin dots */}
          {[370,385,400,415,430].map((x, i) => (
            [100,115,130,145,160].map((y, j) => (
              <circle key={`${i}-${j}`} cx={x} cy={y} r="3" fill="#c4b5fd" opacity="0.7"/>
            ))
          ))}

          <line x1="435" y1="130" x2="475" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr2)"/>

          {/* Suppressor */}
          <rect x="475" y="105" width="80" height="50" rx="6" fill="#fce7f3" stroke="#be185d" strokeWidth="1.5"/>
          <text x="515" y="128" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#be185d">Suppressor</text>
          <text x="515" y="142" textAnchor="middle" fontSize="9" fill="#be185d">↓ background</text>
          <text x="515" y="155" textAnchor="middle" fontSize="9" fill="#be185d">conductivity</text>

          <line x1="555" y1="130" x2="595" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr2)"/>

          {/* Detector */}
          <rect x="595" y="105" width="90" height="50" rx="6" fill="#dcfce7" stroke="#166534" strokeWidth="1.5"/>
          <text x="640" y="126" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#166534">Conductivity</text>
          <text x="640" y="140" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#166534">Detector</text>
          <text x="640" y="154" textAnchor="middle" fontSize="9" fill="#166534">measures ions</text>

          <line x1="685" y1="130" x2="720" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr2)"/>

          {/* Chromatogram box */}
          <rect x="720" y="70" width="120" height="120" rx="6" fill="#fff7ed" stroke="#c2410c" strokeWidth="1.5"/>
          <text x="780" y="90" textAnchor="middle" fontSize="10" fontWeight="bold" fill="#c2410c">Chromatogram</text>
          {/* mini chromatogram peaks */}
          <polyline points="730,165 745,165 750,140 755,165 765,165 773,125 778,165 790,165 800,110 806,165 820,165 830,165" fill="none" stroke="#c2410c" strokeWidth="1.5"/>
          <text x="749" y="180" textAnchor="middle" fontSize="8" fill="#6b7280">F⁻</text>
          <text x="773" y="180" textAnchor="middle" fontSize="8" fill="#6b7280">Cl⁻</text>
          <text x="800" y="180" textAnchor="middle" fontSize="8" fill="#6b7280">SO₄²⁻</text>

          {/* Step labels */}
          <text x="175" y="270" textAnchor="middle" fontSize="10" fill="#6b7280">① Push eluent</text>
          <text x="175" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">at constant flow</text>
          <text x="280" y="270" textAnchor="middle" fontSize="10" fill="#6b7280">② Inject</text>
          <text x="280" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">sample</text>
          <text x="395" y="270" textAnchor="middle" fontSize="10" fill="#6b7280">③ Ions separate</text>
          <text x="395" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">by charge/size</text>
          <text x="515" y="270" textAnchor="middle" fontSize="10" fill="#6b7280">④ Reduce</text>
          <text x="515" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">background</text>
          <text x="640" y="270" textAnchor="middle" fontSize="10" fill="#6b7280">⑤ Detect &amp;</text>
          <text x="640" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">quantify</text>
        </svg>
        <figcaption className="text-center text-sm text-gray-600 mt-2">
          IC schematic for anion analysis: eluent pump → sample injector → ion-exchange column → suppressor → conductivity detector → chromatogram
        </figcaption>
      </figure>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">How it works, step by step</h3>
      <ol className="list-decimal pl-6 mb-6 text-gray-700 space-y-3">
        <li>
          <strong>Eluent and pump.</strong> A high-pressure pump pushes the eluent (mobile phase) through the system at a constant flow rate (typically 1 mL/min). For anion analysis the eluent is typically a dilute sodium carbonate/bicarbonate or sodium hydroxide solution. For cation analysis, dilute sulfuric or methanesulfonic acid is used.
        </li>
        <li>
          <strong>Sample injection.</strong> A measured volume of sample (often 25 μL) is injected via a loop injector. The eluent sweeps it onto the column.
        </li>
        <li>
          <strong>Separation on the column.</strong> The analytical column is packed with polymer beads functionalised with charged groups — quaternary ammonium groups for anion exchange, sulphonate groups for cation exchange. Ions in the sample bind transiently to these sites and are continuously displaced by the eluent ions. Weakly retained ions (like fluoride) elute first; strongly retained ions (like sulphate or phosphate) elute later. Retention time is reproducible and used for identification.
        </li>
        <li>
          <strong>The suppressor — IC&apos;s clever trick.</strong> The conductivity detector measures how well a solution conducts electricity, which depends on the total ion concentration. The problem is that the eluent itself is ionic and would swamp the signal from the analyte ions. The suppressor converts the eluent to water electrochemically: for carbonate eluents it converts Na₂CO₃ into H₂CO₃ (carbonic acid), which barely conducts. Now the detector sees mainly the target analyte ions — greatly boosting sensitivity.
        </li>
        <li>
          <strong>Conductivity detection and quantification.</strong> As each ion band passes the detector, a peak appears on the chromatogram. Peak area (or height) is proportional to concentration. Comparison with external calibration standards gives the concentration of each ion.
        </li>
      </ol>

      <div className="bg-gray-50 border border-gray-200 p-5 my-6">
        <p className="text-gray-700 font-semibold mb-1">Key things to remember about IC</p>
        <ul className="text-gray-700 space-y-1 text-sm">
          <li>— Samples should be aqueous and relatively clean — particulates are filtered (0.22 μm) before injection</li>
          <li>— Very high salt/matrix samples may need dilution to avoid column overloading</li>
          <li>— Detection limits are typically in the μg/L (ppb) to mg/L (ppm) range depending on the analyte</li>
          <li>— Anions and cations require separate columns and eluents — you can&apos;t run both in one injection</li>
          <li>— Common interferences: co-eluting ions, very high chloride suppressing nearby peaks</li>
          <li>— Guard columns protect the analytical column from particulates and strongly retained compounds</li>
        </ul>
      </div>

      {/* ===================== GC ===================== */}
      <h2 className="text-3xl font-bold text-gray-900 mt-14 mb-4 font-mono">GC — Gas Chromatography</h2>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">What is it for?</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        GC separates and quantifies <strong>volatile and semi-volatile organic compounds</strong> — anything that can be vaporised without decomposing. Common applications include solvents and VOCs in air or water, petroleum hydrocarbons, pesticides, flavour and fragrance compounds, and purity testing of organic chemicals. The key requirement is volatility: compounds must exist as a gas at the temperatures used in the column oven (typically up to 350°C).
      </p>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">The core idea: vaporise the sample, then let compounds race through a long thin column</h3>
      <p className="text-gray-700 leading-relaxed mb-4">
        GC uses the same fundamental principle as IC — separation based on differential interaction with a stationary phase — but the mobile phase is an inert carrier gas (helium or nitrogen) rather than a liquid. Inside a very long, narrow capillary column (typically 30–60 m long, 0.25 mm internal diameter) coated on the inside with a thin liquid or polymer stationary phase, compounds partition repeatedly between the gas phase and the stationary phase as they travel through. More volatile compounds spend more time in the gas phase and elute first; heavier, less volatile compounds are retained longer.
      </p>

      {/* GC Diagram */}
      <figure className="my-8">
        <svg viewBox="0 0 840 320" className="w-full max-w-3xl mx-auto" aria-label="Gas chromatography instrument diagram">
          <rect width="840" height="320" fill="#f9fafb" rx="8"/>

          <defs>
            <marker id="arr3" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
              <path d="M0,0 L0,6 L8,3 z" fill="#374151"/>
            </marker>
          </defs>

          {/* Carrier gas */}
          <rect x="15" y="100" width="80" height="60" rx="6" fill="#dbeafe" stroke="#1d4ed8" strokeWidth="1.5"/>
          <text x="55" y="126" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#1e40af">Carrier</text>
          <text x="55" y="140" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#1e40af">Gas</text>
          <text x="55" y="154" textAnchor="middle" fontSize="9" fill="#1e40af">(He or N₂)</text>
          <line x1="95" y1="130" x2="135" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr3)"/>

          {/* Injector */}
          <rect x="135" y="80" width="75" height="100" rx="6" fill="#fef9c3" stroke="#ca8a04" strokeWidth="2"/>
          <text x="172" y="120" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#854d0e">Injector</text>
          <text x="172" y="134" textAnchor="middle" fontSize="9" fill="#854d0e">~250°C</text>
          <text x="172" y="147" textAnchor="middle" fontSize="9" fill="#854d0e">vaporises</text>
          <text x="172" y="160" textAnchor="middle" fontSize="9" fill="#854d0e">sample</text>
          {/* syringe */}
          <line x1="172" y1="80" x2="172" y2="55" stroke="#374151" strokeWidth="2"/>
          <rect x="162" y="42" width="20" height="14" rx="2" fill="#9ca3af" stroke="#374151" strokeWidth="1"/>
          <line x1="172" y1="42" x2="172" y2="30" stroke="#374151" strokeWidth="2"/>
          <text x="195" y="48" textAnchor="start" fontSize="9" fill="#374151">syringe / autosampler</text>

          <line x1="210" y1="130" x2="260" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr3)"/>
          <text x="235" y="122" textAnchor="middle" fontSize="9" fill="#374151">vapour</text>

          {/* Oven with column */}
          <rect x="260" y="60" width="280" height="160" rx="8" fill="#fff7ed" stroke="#c2410c" strokeWidth="2"/>
          <text x="400" y="83" textAnchor="middle" fontSize="12" fontWeight="bold" fill="#c2410c">Column Oven</text>
          <text x="400" y="98" textAnchor="middle" fontSize="10" fill="#c2410c">temperature programmed (e.g. 50°C → 300°C)</text>
          {/* capillary column as coil */}
          <ellipse cx="360" cy="155" rx="55" ry="45" fill="none" stroke="#7e22ce" strokeWidth="3"/>
          <ellipse cx="360" cy="155" rx="40" ry="30" fill="none" stroke="#7e22ce" strokeWidth="3"/>
          <ellipse cx="360" cy="155" rx="25" ry="16" fill="none" stroke="#7e22ce" strokeWidth="3"/>
          <text x="360" y="158" textAnchor="middle" fontSize="9" fill="#7e22ce">capillary</text>
          <text x="360" y="170" textAnchor="middle" fontSize="9" fill="#7e22ce">column</text>
          <text x="360" y="183" textAnchor="middle" fontSize="9" fill="#7e22ce">30–60 m</text>

          {/* stationary phase legend */}
          <rect x="430" y="120" width="95" height="70" rx="4" fill="#f3e8ff" stroke="#7e22ce" strokeWidth="1"/>
          <text x="477" y="138" textAnchor="middle" fontSize="9" fontWeight="bold" fill="#6b21a8">Stationary phase</text>
          <text x="477" y="151" textAnchor="middle" fontSize="8" fill="#6b21a8">coats inner wall</text>
          <text x="477" y="163" textAnchor="middle" fontSize="8" fill="#6b21a8">nonpolar → separates</text>
          <text x="477" y="174" textAnchor="middle" fontSize="8" fill="#6b21a8">by boiling point</text>
          <text x="477" y="185" textAnchor="middle" fontSize="8" fill="#6b21a8">polar → separates</text>
          <text x="477" y="196" textAnchor="middle" fontSize="8" fill="#6b21a8">by polarity too</text>

          <line x1="540" y1="130" x2="590" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr3)"/>
          <text x="565" y="122" textAnchor="middle" fontSize="9" fill="#374151">separated</text>
          <text x="565" y="134" textAnchor="middle" fontSize="9" fill="#374151">compounds</text>

          {/* Detector */}
          <rect x="590" y="95" width="100" height="70" rx="6" fill="#dcfce7" stroke="#166534" strokeWidth="1.5"/>
          <text x="640" y="122" textAnchor="middle" fontSize="11" fontWeight="bold" fill="#166534">Detector</text>
          <text x="640" y="137" textAnchor="middle" fontSize="9" fill="#166534">FID / MS / ECD</text>
          <text x="640" y="150" textAnchor="middle" fontSize="9" fill="#166534">measures signal</text>
          <text x="640" y="163" textAnchor="middle" fontSize="9" fill="#166534">vs. time</text>

          <line x1="690" y1="130" x2="730" y2="130" stroke="#374151" strokeWidth="1.5" markerEnd="url(#arr3)"/>

          {/* Chromatogram */}
          <rect x="730" y="65" width="95" height="130" rx="6" fill="#f0fdf4" stroke="#166534" strokeWidth="1.5"/>
          <text x="777" y="85" textAnchor="middle" fontSize="10" fontWeight="bold" fill="#166534">Output</text>
          <polyline points="738,170 748,170 752,148 756,170 762,170 768,133 773,170 782,170 790,118 796,170 808,170 818,170" fill="none" stroke="#166534" strokeWidth="1.5"/>
          <text x="752" y="183" textAnchor="middle" fontSize="8" fill="#6b7280">A</text>
          <text x="768" y="183" textAnchor="middle" fontSize="8" fill="#6b7280">B</text>
          <text x="791" y="183" textAnchor="middle" fontSize="8" fill="#6b7280">C</text>

          {/* Step labels */}
          <text x="55" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">① Inert</text>
          <text x="55" y="294" textAnchor="middle" fontSize="10" fill="#6b7280">carrier gas</text>
          <text x="172" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">② Vaporise</text>
          <text x="172" y="294" textAnchor="middle" fontSize="10" fill="#6b7280">sample</text>
          <text x="400" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">③ Separate in heated oven</text>
          <text x="400" y="294" textAnchor="middle" fontSize="10" fill="#6b7280">(temperature programmed)</text>
          <text x="640" y="282" textAnchor="middle" fontSize="10" fill="#6b7280">④ Detect &amp;</text>
          <text x="640" y="294" textAnchor="middle" fontSize="10" fill="#6b7280">quantify</text>
        </svg>
        <figcaption className="text-center text-sm text-gray-600 mt-2">
          GC schematic: carrier gas → heated injector → capillary column in programmed oven → detector → chromatogram
        </figcaption>
      </figure>

      <h3 className="text-xl font-bold text-gray-900 mt-6 mb-3">How it works, step by step</h3>
      <ol className="list-decimal pl-6 mb-6 text-gray-700 space-y-3">
        <li>
          <strong>Carrier gas.</strong> An inert gas — usually helium for its low viscosity and good diffusivity, or nitrogen for lower cost — flows continuously through the system at a controlled rate. It carries the vaporised sample through the column but does not interact with it.
        </li>
        <li>
          <strong>Injection and vaporisation.</strong> A small volume of sample (typically 1 μL) is injected via syringe or autosampler into a heated injector port (200–300°C). The sample instantly vaporises and is swept onto the column by the carrier gas. The most common injection mode is <em>split/splitless</em> — in split mode the sample is diluted with gas before entering the column (for concentrated samples); in splitless mode, most of the sample enters the column (for trace analysis).
        </li>
        <li>
          <strong>Separation in the column oven.</strong> The capillary column sits inside an oven whose temperature is precisely programmed — often starting low (to resolve early-eluting volatiles) and ramping up (to push out later compounds in a reasonable time). Inside the column, each compound partitions between the carrier gas and the stationary phase coating the column wall. This partition coefficient determines how long each compound is retained. Compounds elute one at a time, producing discrete peaks.
        </li>
        <li>
          <strong>Detection.</strong> Several detector types are used depending on the application:
          <ul className="list-disc pl-6 mt-2 space-y-1">
            <li><strong>FID (Flame Ionisation Detector)</strong> — the most common. The eluting compounds are burned in a hydrogen/air flame, producing ions that create a measurable current. Responds to almost all organic compounds. Very sensitive and has a wide linear range.</li>
            <li><strong>MS (Mass Spectrometer)</strong> — GC-MS is the gold standard. Each compound is fragmented and its mass spectrum compared to a library — you get positive identification, not just a retention time. Essential for unknown identification and confirmation.</li>
            <li><strong>ECD (Electron Capture Detector)</strong> — extremely sensitive to compounds containing halogens (chlorinated pesticides, PCBs) or nitro groups. Used in environmental trace analysis.</li>
          </ul>
        </li>
        <li>
          <strong>Quantification.</strong> Peak area is proportional to the amount of compound injected. An internal standard (a compound not present in the sample, added at a known concentration) corrects for injection volume variability and instrument drift.
        </li>
      </ol>

      <div className="bg-gray-50 border border-gray-200 p-5 my-6">
        <p className="text-gray-700 font-semibold mb-1">Key things to remember about GC</p>
        <ul className="text-gray-700 space-y-1 text-sm">
          <li>— Compounds must be volatile and thermally stable — non-volatile compounds require derivatisation</li>
          <li>— Water is generally incompatible with FID (damages the flame) — samples are often extracted into an organic solvent first</li>
          <li>— GC-MS allows library matching for positive identification; FID alone only gives a retention time</li>
          <li>— Column choice (polarity of stationary phase) is critical — non-polar columns (e.g. 5% phenyl polysiloxane) separate by boiling point; polar columns separate by polarity too</li>
          <li>— Temperature programming is key to resolving a wide boiling-point range in one run</li>
          <li>— Headspace GC is a common technique: the gas above a heated sample is sampled directly, avoiding matrix effects</li>
        </ul>
      </div>

      {/* ===================== COMPARISON TABLE ===================== */}
      <h2 className="text-3xl font-bold text-gray-900 mt-14 mb-6 font-mono">How They Compare</h2>

      <div className="overflow-x-auto my-6">
        <table className="w-full text-sm border-collapse">
          <thead>
            <tr className="bg-gray-900 text-white">
              <th className="px-4 py-3 text-left font-mono font-bold"></th>
              <th className="px-4 py-3 text-left font-mono font-bold">ICP-OES / ICP-MS</th>
              <th className="px-4 py-3 text-left font-mono font-bold">IC</th>
              <th className="px-4 py-3 text-left font-mono font-bold">GC</th>
            </tr>
          </thead>
          <tbody>
            {[
              ['What it analyses', 'Elements (metals, metalloids)', 'Dissolved ions (anions/cations)', 'Volatile organic compounds'],
              ['Sample state', 'Liquid (aqueous, acidified)', 'Aqueous solution', 'Liquid or gas (must volatilise)'],
              ['Mobile phase', 'Argon gas (carrier)', 'Aqueous eluent (ionic)', 'Inert gas (He or N₂)'],
              ['Separation principle', 'No separation — all elements detected simultaneously', 'Ion-exchange affinity', 'Vapour pressure + stationary phase interaction'],
              ['Common detector', 'CCD (OES) or quadrupole MS', 'Conductivity (suppressed)', 'FID, MS, ECD'],
              ['Detection limits', 'ppb–ppt (OES–MS)', 'ppb–ppm', 'ppb–ppt (with MS/ECD)'],
              ['Gives molecular info?', 'No — elements only', 'Yes — which ion is present', 'Yes — molecular identity (esp. with MS)'],
              ['Typical analytes', 'Pb, Cd, As, Cu, Fe, Ca, P, K…', 'F⁻, Cl⁻, NO₃⁻, SO₄²⁻, PO₄³⁻, Na⁺, K⁺, Ca²⁺', 'Solvents, VOCs, pesticides, hydrocarbons'],
            ].map(([label, icp, ic, gc], i) => (
              <tr key={i} className={i % 2 === 0 ? 'bg-white' : 'bg-gray-50'}>
                <td className="px-4 py-3 font-semibold text-gray-800 border border-gray-200">{label}</td>
                <td className="px-4 py-3 text-gray-700 border border-gray-200">{icp}</td>
                <td className="px-4 py-3 text-gray-700 border border-gray-200">{ic}</td>
                <td className="px-4 py-3 text-gray-700 border border-gray-200">{gc}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* ===================== QC CONTEXT ===================== */}
      <h2 className="text-3xl font-bold text-gray-900 mt-14 mb-4 font-mono">In a QC Lab Context</h2>
      <p className="text-gray-700 leading-relaxed mb-4">
        In a QC environment these techniques sit within a framework of quality assurance. A few things you&apos;ll encounter in practice:
      </p>
      <ul className="list-disc pl-6 mb-6 text-gray-700 space-y-2">
        <li><strong>Calibration curves</strong> — at least three (ideally five or more) calibration standards spanning the expected concentration range, prepared fresh each run. The response should be linear; R² values of ≥0.999 are typically required.</li>
        <li><strong>Blanks</strong> — a reagent blank (solvent/eluent only) checks for contamination; a method blank processes the same reagents through the full sample preparation to catch background contributions.</li>
        <li><strong>QC check standards</strong> — a standard at a known concentration (often a certified reference material) run periodically to confirm the calibration is still valid. Recovery should be within ±10% (or ±20% for trace methods).</li>
        <li><strong>Matrix spikes</strong> — a known amount of analyte added to the sample matrix checks whether the matrix is suppressing or enhancing the signal. Recovery outside 70–130% is a flag for matrix interference.</li>
        <li><strong>Duplicates</strong> — running the same sample twice checks precision. Relative percent difference (RPD) should typically be &lt;20%.</li>
        <li><strong>Retention time windows</strong> (IC, GC) — compound identification relies on retention times matching the standard within a defined window (e.g. ±0.1 min). Drift outside this window requires investigation.</li>
      </ul>

      <h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">Frequently Asked Questions</h2>

      <div className="space-y-6 mb-12">
        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">Why does ICP use argon specifically?</h3>
          <p className="text-gray-700">
            Argon is chemically inert (it won&apos;t react with the sample or torch), it ionises relatively easily to sustain the plasma, and it has a high enough ionisation energy that most analyte elements ionise more readily than the argon itself — keeping the plasma stable. Helium plasmas exist but are much less common.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">What&apos;s the difference between ICP-OES and ICP-MS?</h3>
          <p className="text-gray-700">
            Both use the same plasma source. ICP-OES measures light emitted by excited atoms at characteristic wavelengths — it&apos;s robust, handles high-concentration matrices well, and can run many elements simultaneously. ICP-MS extracts ions from the plasma into a mass spectrometer for separation by mass-to-charge ratio. It achieves detection limits typically 100–1000× lower than OES, and provides isotopic information, but is more sensitive to matrix interferences and more expensive to run.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">Why does IC need a suppressor?</h3>
          <p className="text-gray-700">
            The conductivity detector measures the total conductance of the eluent. Without a suppressor, the ionic eluent (e.g. sodium carbonate) would produce a large background signal that drowns out the analyte peaks. The suppressor converts the eluent counter-ion to a poorly conducting form (e.g. Na⁺ to H⁺ in carbonate, which then becomes water via H⁺ + HCO₃⁻ → H₂CO₃ → H₂O + CO₂). This massively reduces background conductivity and boosts the signal-to-noise ratio for the analyte ions.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">Can GC analyse water samples directly?</h3>
          <p className="text-gray-700">
            Generally no — water is usually removed before GC analysis because it is incompatible with FID detectors and can damage some columns. Options include liquid-liquid extraction (extracting VOCs into an organic solvent), solid-phase extraction, or headspace GC (where the vapour above a heated water sample is sampled directly). Purge-and-trap is another technique that concentrates volatiles from water onto a sorbent before thermal desorption into the GC.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">If I wanted to measure phosphate in a sample, which technique would I use?</h3>
          <p className="text-gray-700">
            It depends on what information you need. IC measures the phosphate ion (PO₄³⁻) directly and gives you the inorganic phosphate concentration. ICP-OES or ICP-MS measures total phosphorus — the sum of all forms of phosphorus in the sample (inorganic phosphate, organic phosphate compounds, polyphosphates, etc.) after acid digestion. For fertiliser or water QC where you need to know the ionic species specifically, IC is often preferred; for total elemental content, ICP is the choice.
          </p>
        </div>

        <div className="border-l-4 border-gray-400 pl-6">
          <h3 className="text-xl font-bold text-gray-900 mb-2">What does &quot;matrix effect&quot; mean in analytical chemistry?</h3>
          <p className="text-gray-700">
            The matrix is everything in a sample apart from the target analyte — salts, organic matter, acids, other ions. Matrix effects occur when these background constituents interfere with the measurement: they might suppress the signal (e.g. high dissolved solids reducing nebulisation efficiency in ICP), enhance it (e.g. organic solvents increasing sensitivity in ICP-MS), or co-elute with analytes in IC or GC. Managing matrix effects is a large part of method development — solutions include dilution, matrix-matched calibration, internal standards, and sample clean-up steps.
          </p>
        </div>
      </div>

      {/* Author byline */}
      <div className="mt-12 pt-8 border-t border-gray-200">
        <p className="text-gray-600 italic">
          Written by <strong>Miles Christou</strong>, chemistry graduate with practical experience in analytical technique development including HPLC-MS and environmental sample preparation
        </p>
      </div>
    </div>
  );
}
