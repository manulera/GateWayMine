import { Box, CircularProgress, Typography } from '@mui/material';
import React from 'react'
import { useState } from 'react'
import FilterBar from './FilterBar';
import ResultsTable from './ResultsTable';
import PlasmidTableRowContent from './PlasmidTableRowContent';
import './app.css';

function App() {

  const [data, setData] = useState(null);
  const [sites, setSites] = useState([]);
  const [features, setFeatures] = useState([]);
  const [kits, setKits] = useState([]);

  const [selectedSites, setSelectedSites] = useState([]);
  const [selectedFeatures, setSelectedFeatures] = useState([]);
  const [selectedSource, setSelectedSource] = useState('all');
  const [results, setResults] = useState([]);
  const [nameFilter, setNameFilter] = useState('');
  const [kitFilter, setKitFilter] = useState('');


  React.useEffect(() => {
    const updateData = async () => {
      try {
        const response = await fetch('/plasmid_features.json');
        const dataUnfiltered = await response.json();
        const uniqueSites = new Set();
        const uniqueFeatures = new Set();
        const uniqueKits = new Set();
        const dataFiltered = dataUnfiltered.filter((plasmid) => {
          return plasmid.att_sites && plasmid.features
        })
        dataFiltered.forEach((plasmid) => {
          uniqueSites.add(...plasmid.att_sites)
          uniqueFeatures.add(...plasmid.features)
          if (plasmid.kit) {
            uniqueKits.add(plasmid.kit.name)
          }
        })
        setSites([...uniqueSites].sort());
        setFeatures([...uniqueFeatures].sort());
        setKits([...uniqueKits].sort());
        // Pre-render the rows for performance
        dataFiltered.forEach((plasmid, index) => {
          plasmid.row = PlasmidTableRowContent({ row: plasmid });
        })
        // Sort the addgene rows to the top
        dataFiltered.sort((a, b) => a.source === 'addgene' ? -1 : b.source === 'addgene' ? 1 : 0);
        setData(dataFiltered);
      } catch (error) {
        console.error('Error fetching data:', error);
      }
    }
    updateData();
  }, []);

  React.useEffect(() => {
    if (selectedSites.length === 0 && selectedFeatures.length === 0 && nameFilter.length < 2 && kitFilter === '') {
      setResults([]);
      return;
    }
    const filteredResults = data.filter(plasmid => {
      const matchesSites = selectedSites.length === 0 ||
        selectedSites.every(site => plasmid.att_sites.includes(site));

      const matchesFeatures = selectedFeatures.length === 0 ||
        selectedFeatures.every(feature => plasmid.features.includes(feature));

      const query = nameFilter.toLowerCase();
      const matchesName = query === '' ||
        plasmid.plasmid_name.toLowerCase().includes(query);

      const matchesSource = selectedSource === 'all' || plasmid.source === selectedSource;

      const matchesKit = kitFilter === '' || plasmid.kit?.name === kitFilter;

      return matchesSites && matchesFeatures && matchesName && matchesSource && matchesKit;
    });
    const filteredFeatures = new Set(filteredResults.flatMap(plasmid => plasmid.features));
    const filteredSites = new Set(filteredResults.flatMap(plasmid => plasmid.att_sites));
    console.log(filteredResults);
    console.log(filteredFeatures);
    setResults(filteredResults);
    // Remove selected features that are not in the filtered features
    const newSelectedFeatures = selectedFeatures.filter(feature => filteredFeatures.has(feature));
    if (newSelectedFeatures.length !== selectedFeatures.length) {
      setSelectedFeatures(newSelectedFeatures);
    }
    // Remove selected sites that are not in the filtered sites 
    const newSelectedSites = selectedSites.filter(site => filteredSites.has(site));
    if (newSelectedSites.length !== selectedSites.length) {
      setSelectedSites(newSelectedSites);
    }

    setFeatures([...filteredFeatures]);
    setSites([...filteredSites]);
  }, [selectedSites, selectedFeatures, nameFilter, selectedSource, kitFilter])

  const title = (
    <Typography gutterBottom textAlign="center">
      <h1 style={{ fontSize: '4em', fontWeight: 'normal', marginBottom: 0, paddingBottom: 0 }}>GateWayMine</h1>
      <h2 style={{ fontWeight: 'normal', marginTop: 0 }}>A database of Gateway plasmids.</h2>
      <p>Visit the <a href="https://github.com/manulera/gateway_sequences" target="_blank" rel="noopener noreferrer">GitHub repository</a></p>
    </Typography>
  )

  if (!data) {
    return (
      <>
        {title}
        <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '50vh' }}>
          <CircularProgress />
        </Box>
      </>
    )
  }

  return (
    <>
      {title}
      <Typography gutterBottom textAlign="center">
        <Box sx={{ p: 2, display: 'flex', justifyContent: 'center', alignItems: 'center', fontSize: '1.2rem' }}>
          <FilterBar
            sites={sites}
            features={features}
            setSelectedSites={setSelectedSites}
            setSelectedFeatures={setSelectedFeatures}
            nameFilter={nameFilter}
            setNameFilter={setNameFilter}
            selectedSource={selectedSource}
            setSelectedSource={setSelectedSource}
            kits={kits}
            kitFilter={kitFilter}
            setKitFilter={setKitFilter}
          />
        </Box>
        {results.length > 0 && (
          <div>
            {results.length} results
          </div>
        )}
        <ResultsTable results={results} />
      </Typography>
    </>
  )
}

export default App
