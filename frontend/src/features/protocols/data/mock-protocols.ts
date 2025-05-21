/**
 * Mock data for protocols
 */
export const MOCK_PROTOCOLS = [
  {
    id: '1',
    name: 'Standard Cell Freezing Protocol',
    version: '2.1.0',
    description: 'Basic protocol for freezing mammalian cells using DMSO as a cryoprotectant with controlled rate cooling.',
    author: { name: 'Dr. John Doe', id: 'user1' },
    created_at: '2023-09-01',
    updated_at: '2023-10-15',
    duration: '2 hours',
    is_template: false,
    cell_types: ['Human embryonic stem cells', 'Primary fibroblasts', 'CHO cells'],
    equipment: [
      'Controlled-rate freezer',
      'Cryovials',
      'Liquid nitrogen storage tank',
      'Water bath (37°C)',
      'Cell culture materials'
    ],
    cryoprotectants: [
      { id: '1', name: 'DMSO', concentration: '10%' },
      { id: '2', name: 'Fetal bovine serum', concentration: '20%' }
    ],
    freezing_rate: '-1°C/min',
    storage_temperature: '-196°C (liquid nitrogen)',
    thawing_method: 'Rapid thawing in 37°C water bath',
    steps: [
      {
        id: 's1',
        order: 1,
        title: 'Prepare cryopreservation medium',
        description: 'Mix cell culture medium with DMSO to achieve 10% final concentration and FBS to 20% final concentration. Keep on ice.',
        duration: '15 minutes',
        is_critical: true
      },
      {
        id: 's2',
        order: 2,
        title: 'Prepare cells',
        description: 'Harvest cells in log phase growth (70-80% confluency). Count cells and check viability.',
        duration: '20 minutes',
        is_critical: true
      },
      {
        id: 's3',
        order: 3,
        title: 'Centrifuge cell suspension',
        description: 'Centrifuge at 200 x g for 5 minutes at room temperature.',
        duration: '10 minutes',
        is_critical: false
      },
      {
        id: 's4',
        order: 4,
        title: 'Resuspend cells in cryopreservation medium',
        description: 'Gently resuspend cell pellet in cold cryopreservation medium to achieve 1-5 x 10^6 cells/mL.',
        duration: '10 minutes',
        is_critical: true
      },
      {
        id: 's5',
        order: 5,
        title: 'Aliquot into cryovials',
        description: 'Transfer 1 mL of cell suspension to each pre-labeled cryovial.',
        duration: '15 minutes',
        is_critical: false
      },
      {
        id: 's6',
        order: 6,
        title: 'Controlled-rate freezing',
        description: 'Place cryovials in controlled-rate freezer and run program: hold at 4°C for 10 minutes, then cool at -1°C/min to -80°C.',
        duration: '90 minutes',
        is_critical: true
      },
      {
        id: 's7',
        order: 7,
        title: 'Transfer to long-term storage',
        description: 'Quickly transfer frozen vials to liquid nitrogen storage tank.',
        duration: '10 minutes',
        is_critical: true
      }
    ],
    notes: 'For optimal results, cells should be in log phase of growth and have viability >90% before freezing. DMSO is toxic to cells at room temperature, so work quickly when adding cells to freezing medium.',
    references: [
      { title: 'Optimization of cryopreservation procedures for mammalian cells', authors: 'Smith et al.', journal: 'Cryobiology', year: 2020 },
      { title: 'Best practices for cell banking', authors: 'Johnson & Williams', journal: 'Nature Methods', year: 2019 }
    ],
    used_in_experiments: [
      { id: '1', title: 'DMSO concentration optimization', date: '2023-11-10' },
      { id: '3', title: 'Fibroblast freezing comparison', date: '2023-10-20' }
    ]
  },
  {
    id: '2',
    name: 'Vitrification Protocol',
    version: '2.0.1',
    description: 'Advanced vitrification method with controlled cooling rates and combination cryoprotectant solution.',
    author: { name: 'Dr. Jane Smith', id: 'user2' },
    created_at: '2023-08-15',
    updated_at: '2023-10-20',
    duration: '3.5 hours',
    is_template: false,
    cell_types: ['Embryos', 'Oocytes', 'Stem cells'],
    equipment: [
      'Open pulled straws',
      'Cryotop device',
      'Liquid nitrogen container',
      'Water bath (37°C)',
      'Stereomicroscope'
    ],
    cryoprotectants: [
      { id: '4', name: 'Propylene glycol', concentration: '15%' },
      { id: '5', name: 'Ethylene glycol', concentration: '15%' },
      { id: '3', name: 'Trehalose', concentration: '0.5M' }
    ],
    freezing_rate: 'Ultra-rapid (>20,000°C/min)',
    storage_temperature: '-196°C (liquid nitrogen)',
    thawing_method: 'Ultra-rapid warming in 37°C medium with stepwise dilution',
    steps: [
      {
        id: 's1',
        order: 1,
        title: 'Prepare equilibration solution',
        description: 'Mix base medium with 7.5% ethylene glycol and 7.5% propylene glycol.',
        duration: '20 minutes',
        is_critical: true
      },
      {
        id: 's2',
        order: 2,
        title: 'Prepare vitrification solution',
        description: 'Mix base medium with 15% ethylene glycol, 15% propylene glycol, and 0.5M trehalose.',
        duration: '20 minutes',
        is_critical: true
      },
      {
        id: 's3',
        order: 3,
        title: 'Equilibration step',
        description: 'Expose cells/tissues to equilibration solution for 10 minutes at room temperature.',
        duration: '15 minutes',
        is_critical: true
      },
      {
        id: 's4',
        order: 4,
        title: 'Vitrification step',
        description: 'Transfer to vitrification solution for 1 minute at room temperature.',
        duration: '5 minutes',
        is_critical: true
      },
      {
        id: 's5',
        order: 5,
        title: 'Load onto carrier',
        description: 'Place cells/tissues onto Cryotop device with minimal volume of vitrification solution.',
        duration: '5 minutes',
        is_critical: true
      },
      {
        id: 's6',
        order: 6,
        title: 'Plunge into liquid nitrogen',
        description: 'Rapidly immerse loaded Cryotop into liquid nitrogen.',
        duration: '5 minutes',
        is_critical: true
      },
      {
        id: 's7',
        order: 7,
        title: 'Transfer to storage',
        description: 'Place Cryotop into storage container filled with liquid nitrogen.',
        duration: '10 minutes',
        is_critical: true
      }
    ],
    notes: 'This protocol is specifically optimized for small volume samples. Timing is critical, especially in the vitrification step to prevent cytotoxicity while ensuring sufficient cryoprotectant permeation.',
    references: [
      { title: 'Advances in vitrification methods for cryopreservation of embryos and oocytes', authors: 'Kuwayama et al.', journal: 'Reproductive BioMedicine Online', year: 2022 },
      { title: 'Optimization of cryoprotectant solutions for vitrification', authors: 'Chen & Li', journal: 'Cryobiology', year: 2021 }
    ],
    used_in_experiments: [
      { id: '5', title: 'Ethylene Glycol Protocol', date: '2023-11-05' }
    ]
  },
  {
    id: '3',
    name: 'Plant Tissue Preservation',
    version: '1.5.0',
    description: 'Specialized protocol for preserving plant tissue samples using glycerol as the primary cryoprotectant.',
    author: { name: 'Dr. Sarah Johnson', id: 'user3' },
    created_at: '2023-07-10',
    updated_at: '2023-09-25',
    duration: '4 hours',
    is_template: true,
    cell_types: ['Arabidopsis leaf cells', 'Plant meristems', 'Callus cultures'],
    equipment: [
      'Programmable freezer',
      'Cryovials',
      'Ultra-low freezer (-80°C)',
      'Sterile work cabinet',
      'Plant growth chamber'
    ],
    cryoprotectants: [
      { id: '2', name: 'Glycerol', concentration: '15%' },
      { id: '3', name: 'Trehalose', concentration: '0.4M' }
    ],
    freezing_rate: '-0.5°C/min',
    storage_temperature: '-80°C (ultra-low freezer)',
    thawing_method: 'Controlled rate thawing with stepwise dilution',
    steps: [
      {
        id: 's1',
        order: 1,
        title: 'Prepare cryoprotectant solution',
        description: 'Mix plant culture medium with 15% glycerol and 0.4M trehalose. Filter sterilize.',
        duration: '30 minutes',
        is_critical: true
      },
      {
        id: 's2',
        order: 2,
        title: 'Prepare plant tissues',
        description: 'Harvest actively growing tissue samples of approximately 2-3mm in size.',
        duration: '45 minutes',
        is_critical: true
      },
      {
        id: 's3',
        order: 3,
        title: 'Preconditioning',
        description: 'Place tissues on osmotic medium (0.3M sucrose) for 24h at 4°C.',
        duration: '24 hours',
        is_critical: true
      },
      {
        id: 's4',
        order: 4,
        title: 'Cryoprotectant loading',
        description: 'Transfer tissues to cryoprotectant solution for 60 minutes at 4°C with gentle agitation.',
        duration: '60 minutes',
        is_critical: true
      },
      {
        id: 's5',
        order: 5,
        title: 'Transfer to cryovials',
        description: 'Place tissues in cryovials with 0.5ml of cryoprotectant solution.',
        duration: '20 minutes',
        is_critical: false
      },
      {
        id: 's6',
        order: 6,
        title: 'Controlled-rate freezing',
        description: 'Place cryovials in programmable freezer. Cool at -0.5°C/min to -40°C, hold for 30min, then cool at -5°C/min to -80°C.',
        duration: '120 minutes',
        is_critical: true
      },
      {
        id: 's7',
        order: 7,
        title: 'Transfer to storage',
        description: 'Transfer cryovials to -80°C ultra-low freezer for long-term storage.',
        duration: '10 minutes',
        is_critical: true
      }
    ],
    notes: 'This protocol is optimized for preservation of plant genetic resources. Preconditioning is essential for successful cryopreservation of most plant tissues.',
    references: [
      { title: 'Cryopreservation of plant genetic resources', authors: 'Panis & Lambardi', journal: 'Plant Cell Reports', year: 2018 },
      { title: 'Advances in plant cryopreservation', authors: 'Reed et al.', journal: 'CryoLetters', year: 2020 }
    ],
    used_in_experiments: [
      { id: '3', title: 'Glycerol Cryopreservation', date: '2023-10-15' }
    ]
  }
];