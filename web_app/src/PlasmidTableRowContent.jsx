import React from 'react'
import { TableRow, TableCell, Box, Chip } from '@mui/material'

function PlasmidLink(plasmid) {
    if (plasmid.source === 'snapgene') {
        return "https://www.snapgene.com/plasmids"
    }
    if (plasmid.source === 'addgene') {
        return `https://www.addgene.org/${plasmid.addgene_id}`
    }
}

function formatReference(reference) {
    if (reference.includes("doi.org")) {
        return `doi:${reference.split("doi.org/")[1]}`
    }
    if (reference.includes("pubmed")) {
        return `PubMed:${reference.split("pubmed/")[1]}`
    }
    throw new Error(`Unknown reference format: ${reference}`)

}

function PlasmidTableRowContent({ row }) {

    return (
        <>
            <TableCell>{row.plasmid_name}</TableCell>
            <TableCell>
                <Chip
                    label={row.source}
                    component="a"
                    target="_blank"
                    rel="noopener noreferrer"
                    href={PlasmidLink(row)}
                    clickable
                    size="small"
                    sx={{ mt: 0.5 }}
                    color="primary"
                />
                {row.references?.map((ref, i) => (
                    <Chip
                        label={formatReference(ref)}
                        component="a"
                        target="_blank"
                        rel="noopener noreferrer"
                        href={ref}
                        clickable
                        size="small"
                        sx={{ mt: 0.5 }}
                        color="success"
                    />
                ))}
                {row.kit && <Chip
                    label="kit"
                    component="a"
                    target="_blank"
                    rel="noopener noreferrer"
                    href={row.kit.url}
                    clickable
                    size="small"
                    sx={{ mt: 0.5 }}
                    color="secondary"
                />}
            </TableCell>
            <TableCell>
                <Box sx={{ display: 'flex', gap: 0.5, flexWrap: 'wrap' }}>
                    {row.att_sites.map((site, i) => (
                        <Chip key={i} label={site} size="small" />
                    ))}
                </Box>
            </TableCell>
            <TableCell>
                <Box sx={{ display: 'flex', gap: 0.5, flexWrap: 'wrap' }}>
                    {row.features.map((feature, i) => (
                        <Chip key={i} label={feature} size="small" />
                    ))}
                </Box>
            </TableCell>
        </>
    )
}

export default PlasmidTableRowContent