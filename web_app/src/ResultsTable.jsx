import React from 'react'
import { Table, TableBody, TableCell, TableContainer, TableHead, TableRow, Paper, TablePagination, Chip, Box } from '@mui/material'

function ResultsTable({ results }) {
    const [page, setPage] = React.useState(0);
    const [rowsPerPage] = React.useState(50);

    React.useEffect(() => {
        setPage(0);
    }, [results]);

    const handleChangePage = (event, newPage) => {
        setPage(newPage);
    };

    return (
        <Paper sx={{ width: '100%', maxWidth: 1200, mx: 'auto', mt: 2 }}>
            <TableContainer>
                <Table>
                    <TableHead>
                        <TableRow>
                            <TableCell>Plasmid</TableCell>
                            <TableCell>Links</TableCell>
                            <TableCell>Sites</TableCell>
                            <TableCell>Features</TableCell>
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {results
                            .slice(page * rowsPerPage, page * rowsPerPage + rowsPerPage)
                            .map((row, index) => (
                                <TableRow key={index} sx={{
                                    '&:hover': {
                                        backgroundColor: 'rgba(0, 0, 0, 0.04)',
                                    },
                                    backgroundColor: index % 2 === 0 ? 'inherit' : 'rgba(0, 0, 0, 0.02)'
                                }} >
                                    {row.row}
                                </TableRow>
                            ))}
                    </TableBody>
                </Table>
            </TableContainer>
            <TablePagination
                component="div"
                count={results.length}
                rowsPerPage={rowsPerPage}
                page={page}
                onPageChange={handleChangePage}
                rowsPerPageOptions={[100]}
            />
        </Paper>
    )
}

export default ResultsTable