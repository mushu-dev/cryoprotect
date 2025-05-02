/**
 * CryoProtect Analyzer Web Interface
 * Charts Module - Handles data visualization using Chart.js
 */

const Charts = (function() {
  // Chart color palette
  const COLORS = [
    '#2c3e50', // Primary
    '#3498db', // Secondary
    '#e74c3c', // Accent
    '#2ecc71', // Success
    '#f39c12', // Warning
    '#9b59b6', // Purple
    '#1abc9c', // Teal
    '#34495e', // Dark
    '#e67e22', // Orange
    '#95a5a6'  // Gray
  ];
  
  /**
   * Create a chart for molecule properties
   * @param {string} canvasId - The ID of the canvas element
   * @param {Array} molecules - Array of molecule objects
   */
  function createMoleculePropertiesChart(canvasId, molecules) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) {
      console.error(`Canvas element with ID ${canvasId} not found`);
      return;
    }
    
    // Extract top molecules by total score
    const topMolecules = molecules
      .filter(molecule => molecule.properties && molecule.properties['Total Score'])
      .sort((a, b) => b.properties['Total Score'] - a.properties['Total Score'])
      .slice(0, 10);
    
    // Prepare data for chart
    const labels = topMolecules.map(molecule => molecule.name || `CID: ${molecule.cid}`);
    const scores = topMolecules.map(molecule => molecule.properties['Total Score']);
    
    // Create chart
    new Chart(canvas, {
      type: 'bar',
      data: {
        labels: labels,
        datasets: [{
          label: 'Total Score',
          data: scores,
          backgroundColor: COLORS[0],
          borderColor: COLORS[0],
          borderWidth: 1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        aspectRatio: function() {
          // Adjust aspect ratio based on screen width
          return window.innerWidth < 768 ? 1 : 2;
        },
        onResize: function(chart, size) {
          // Dynamically adjust font sizes on resize
          const fontSize = size.width < 768 ? 14 : size.width < 992 ? 16 : 18;
          chart.options.plugins.title.font.size = fontSize;
          chart.options.scales.x.ticks.font.size = fontSize - 2;
          chart.options.scales.y.ticks.font.size = fontSize - 2;
        },
        plugins: {
          title: {
            display: true,
            text: 'Top Molecules by Total Score',
            font: {
              size: 16,
              weight: 'bold'
            },
            padding: {
              top: 10,
              bottom: 15
            }
          },
          legend: {
            display: false,
            position: function() {
              // Adjust legend position based on screen width
              return window.innerWidth < 768 ? 'bottom' : 'top';
            },
            labels: {
              boxWidth: function() {
                return window.innerWidth < 768 ? 10 : 15;
              },
              padding: function() {
                return window.innerWidth < 768 ? 8 : 15;
              }
            }
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                return `Score: ${context.raw}`;
              }
            },
            // Enhanced touch interaction for mobile
            usePointStyle: true,
            boxPadding: 6,
            bodyFont: {
              size: function() {
                return window.innerWidth < 768 ? 12 : 14;
              }
            }
          }
        },
        scales: {
          y: {
            beginAtZero: true,
            ticks: {
              font: {
                size: 12
              },
              // Ensure enough ticks on small screens but not too many
              maxTicksLimit: function() {
                return window.innerWidth < 768 ? 5 : 8;
              }
            },
            title: {
              display: true,
              text: 'Total Score'
            }
          },
          x: {
            ticks: {
              font: {
                size: 12
              },
              // Dynamically adjust rotation based on screen width and label length
              maxRotation: function() {
                return window.innerWidth < 768 ? 90 : 45;
              },
              minRotation: function() {
                return window.innerWidth < 768 ? 45 : 0;
              },
              autoSkip: true,
              autoSkipPadding: 10,
              maxTicksLimit: function() {
                return window.innerWidth < 768 ? 5 : 10;
              }
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a chart for molecule property distribution
   * @param {string} canvasId - The ID of the canvas element
   * @param {Array} molecules - Array of molecule objects
   * @param {string} propertyName - The name of the property to visualize
   */
  function createPropertyDistributionChart(canvasId, molecules, propertyName) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) {
      console.error(`Canvas element with ID ${canvasId} not found`);
      return;
    }
    
    // Extract property values
    const propertyValues = molecules
      .filter(molecule => molecule.properties && molecule.properties[propertyName] !== undefined)
      .map(molecule => molecule.properties[propertyName]);
    
    if (propertyValues.length === 0) {
      console.error(`No data found for property ${propertyName}`);
      return;
    }
    
    // For numeric properties, create a histogram
    if (typeof propertyValues[0] === 'number') {
      createHistogram(canvas, propertyValues, propertyName);
    } 
    // For text properties, create a pie chart
    else if (typeof propertyValues[0] === 'string') {
      createPieChart(canvas, propertyValues, propertyName);
    }
  }
  
  /**
   * Create a histogram for numeric property values
   * @param {HTMLElement} canvas - The canvas element
   * @param {Array} values - Array of numeric values
   * @param {string} propertyName - The name of the property
   */
  function createHistogram(canvas, values, propertyName) {
    // Calculate histogram bins
    const min = Math.min(...values);
    const max = Math.max(...values);
    const binCount = 10;
    const binWidth = (max - min) / binCount;
    
    const bins = Array(binCount).fill(0);
    values.forEach(value => {
      const binIndex = Math.min(Math.floor((value - min) / binWidth), binCount - 1);
      bins[binIndex]++;
    });
    
    // Create bin labels
    const labels = [];
    for (let i = 0; i < binCount; i++) {
      const binStart = min + i * binWidth;
      const binEnd = min + (i + 1) * binWidth;
      labels.push(`${binStart.toFixed(1)} - ${binEnd.toFixed(1)}`);
    }
    
    // Create chart
    new Chart(canvas, {
      type: 'bar',
      data: {
        labels: labels,
        datasets: [{
          label: 'Frequency',
          data: bins,
          backgroundColor: COLORS[1],
          borderColor: COLORS[1],
          borderWidth: 1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: `Distribution of ${propertyName}`,
            font: {
              size: 16
            }
          },
          legend: {
            display: false
          }
        },
        scales: {
          y: {
            beginAtZero: true,
            title: {
              display: true,
              text: 'Frequency'
            }
          },
          x: {
            title: {
              display: true,
              text: propertyName
            },
            ticks: {
              maxRotation: 45,
              minRotation: 45
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a pie chart for text property values
   * @param {HTMLElement} canvas - The canvas element
   * @param {Array} values - Array of text values
   * @param {string} propertyName - The name of the property
   */
  function createPieChart(canvas, values, propertyName) {
    // Count occurrences of each value
    const valueCounts = {};
    values.forEach(value => {
      valueCounts[value] = (valueCounts[value] || 0) + 1;
    });
    
    // Prepare data for chart
    const labels = Object.keys(valueCounts);
    const data = Object.values(valueCounts);
    
    // Create chart
    new Chart(canvas, {
      type: 'pie',
      data: {
        labels: labels,
        datasets: [{
          data: data,
          backgroundColor: COLORS.slice(0, labels.length),
          borderWidth: 1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: `Distribution of ${propertyName}`,
            font: {
              size: 16
            }
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const value = context.raw;
                const total = context.dataset.data.reduce((a, b) => a + b, 0);
                const percentage = ((value / total) * 100).toFixed(1);
                return `${context.label}: ${value} (${percentage}%)`;
              }
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a chart for mixture composition
   * @param {string} canvasId - The ID of the canvas element
   * @param {Array} mixtures - Array of mixture objects
   */
  function createMixtureCompositionChart(canvasId, mixtures) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) {
      console.error(`Canvas element with ID ${canvasId} not found`);
      return;
    }
    
    // Get the first mixture with components
    const mixture = mixtures.find(mix => mix.components && mix.components.length > 0);
    
    if (!mixture) {
      console.error('No mixture with components found');
      return;
    }
    
    // Prepare data for chart
    const labels = mixture.components.map(component => component.name || `CID: ${component.cid}`);
    const data = mixture.components.map(component => component.concentration);
    
    // Create chart
    new Chart(canvas, {
      type: 'pie',
      data: {
        labels: labels,
        datasets: [{
          data: data,
          backgroundColor: COLORS.slice(0, labels.length),
          borderWidth: 1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: `Composition of ${mixture.name}`,
            font: {
              size: 16
            }
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const value = context.raw;
                const total = context.dataset.data.reduce((a, b) => a + b, 0);
                const percentage = ((value / total) * 100).toFixed(1);
                const unit = mixture.components[context.dataIndex].concentration_unit;
                return `${context.label}: ${value}${unit} (${percentage}%)`;
              }
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a chart comparing prediction with experiment
   * @param {string} canvasId - The ID of the canvas element
   * @param {Object} comparison - Comparison data
   */
  function createComparisonChart(canvasId, comparison) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) {
      console.error(`Canvas element with ID ${canvasId} not found`);
      return;
    }
    
    if (!comparison.prediction || !comparison.experiment) {
      console.error('Invalid comparison data');
      return;
    }
    
    // Prepare data for chart
    const labels = ['Prediction', 'Experiment'];
    const data = [
      comparison.prediction.numeric_value,
      comparison.experiment.numeric_value
    ];
    
    // Create chart
    new Chart(canvas, {
      type: 'bar',
      data: {
        labels: labels,
        datasets: [{
          label: 'Value',
          data: data,
          backgroundColor: [COLORS[1], COLORS[2]],
          borderColor: [COLORS[1], COLORS[2]],
          borderWidth: 1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: 'Prediction vs. Experiment',
            font: {
              size: 16
            }
          },
          tooltip: {
            callbacks: {
              afterLabel: function(context) {
                if (context.dataIndex === 0) {
                  return `Confidence: ${(comparison.prediction.confidence * 100).toFixed(1)}%`;
                } else {
                  return `Conditions: ${comparison.experiment.conditions || 'Not specified'}`;
                }
              }
            }
          }
        },
        scales: {
          y: {
            beginAtZero: true,
            title: {
              display: true,
              text: 'Value'
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a radar chart for molecule properties
   * @param {string} canvasId - The ID of the canvas element
   * @param {Object} molecule - Molecule object
   */
  function createMoleculeRadarChart(canvasId, molecule) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) {
      console.error(`Canvas element with ID ${canvasId} not found`);
      return;
    }
    
    if (!molecule.properties) {
      console.error('Molecule has no properties');
      return;
    }
    
    // Select properties for radar chart
    const propertyNames = [
      'Hydrogen Bonding Score',
      'Solubility Score',
      'Membrane Permeability Score',
      'Toxicity Score',
      'Stability Score',
      'Environmental Score'
    ];
    
    // Filter properties that exist in the molecule
    const availableProperties = propertyNames.filter(name => 
      molecule.properties[name] !== undefined
    );
    
    if (availableProperties.length < 3) {
      console.error('Not enough properties for radar chart');
      return;
    }
    
    // Prepare data for chart
    const data = availableProperties.map(name => molecule.properties[name]);
    
    // Create chart
    new Chart(canvas, {
      type: 'radar',
      data: {
        labels: availableProperties,
        datasets: [{
          label: molecule.name || `CID: ${molecule.cid}`,
          data: data,
          backgroundColor: `${COLORS[0]}40`,
          borderColor: COLORS[0],
          borderWidth: 2,
          pointBackgroundColor: COLORS[0],
          pointRadius: 4
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: 'Property Profile',
            font: {
              size: 16
            }
          }
        },
        scales: {
          r: {
            angleLines: {
              display: true
            },
            suggestedMin: 0,
            suggestedMax: 100
          }
        }
      }
    });
  }
  
  // Public API
  return {
    createMoleculePropertiesChart,
    createPropertyDistributionChart,
    createMixtureCompositionChart,
    createComparisonChart,
    createMoleculeRadarChart
  };
})();