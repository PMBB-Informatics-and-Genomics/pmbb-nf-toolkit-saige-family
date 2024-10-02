document.addEventListener('DOMContentLoaded', () => {
    console.log('Script loaded and DOM fully parsed');

    let selectedCohort = '';
    let selectedPhenotype = '';
    let selectedGroupTest = '';
    let selectedMaf = '';
    let allData = [];
    let filteredData = [];
    let columns = [];
    let currentFilters = {};
    let uniqueValues = {};
    let plotManifest = [];
    let cohorts = [];


    function filterTable() {
        filteredData = allData.filter(row => {
            return columns.every(col => {
                const columnValue = (row[col] || '').toString().toUpperCase();
                return columnValue.includes((currentFilters[col] || '').toUpperCase());
            });
        });
        renderTableBody(filteredData);
    }

    window.toggleSubmenu = function(id) {
        const submenu = document.getElementById(id);
        if (submenu.style.display === 'none' || submenu.style.display === '') {
            submenu.style.display = 'block';
        } else {
            submenu.style.display = 'none';
        }
    }

    window.showAnalysisDescriptions = function(event) {
        event.preventDefault();
        document.getElementById('default-view').style.display = 'none';
        document.getElementById('statistics').style.display = 'none';
        document.querySelector('.table-container').style.display = 'none';
        document.getElementById('plot-tabs').style.display = 'none';
        document.getElementById('cohort-buttons').style.display = 'none';
        document.getElementById('hits-table-title').style.display = 'none';
        document.getElementById('phenotype-filter').style.display = 'none';
        document.getElementById('analysis-descriptions').style.display = 'block';
        toggleSubmenu('analysis-description-submenu');
    }

    window.showStatsDescription = function(event) {
        event.preventDefault();
        document.getElementById('default-view').style.display = 'none';
        document.getElementById('statistics').style.display = 'none';
        document.querySelector('.table-container').style.display = 'none';
        document.getElementById('plot-tabs').style.display = 'none';
        document.getElementById('cohort-buttons').style.display = 'none';
        document.getElementById('hits-table-title').style.display = 'none';
        document.getElementById('phenotype-filter').style.display = 'none';
        document.getElementById('analysis-descriptions').style.display = 'block';
        document.getElementById('stats-description').style.display = 'block';
        document.getElementById('method-summary').style.display = 'none';
    }

    window.showMethodSummary = function(event) {
        event.preventDefault();
        document.getElementById('default-view').style.display = 'none';
        document.getElementById('statistics').style.display = 'none';
        document.querySelector('.table-container').style.display = 'none';
        document.getElementById('plot-tabs').style.display = 'none';
        document.getElementById('cohort-buttons').style.display = 'none';
        document.getElementById('hits-table-title').style.display = 'none';
        document.getElementById('phenotype-filter').style.display = 'none';
        document.getElementById('analysis-descriptions').style.display = 'block';
        document.getElementById('stats-description').style.display = 'none';
        document.getElementById('method-summary').style.display = 'block';
    }

    window.showCohort = function(cohort) {
        selectedCohort = cohort;
        document.getElementById('default-view').style.display = 'none';
        document.getElementById('statistics').style.display = 'none';
        document.querySelector('.table-container').style.display = 'block';
        document.getElementById('plot-tabs').style.display = 'block';
        document.getElementById('plot-container').style.display = 'block';
        document.getElementById('cohort-buttons').style.display = 'none';
        document.getElementById('hits-table-title').style.display = 'block';
        document.getElementById('phenotype-filter').style.display = 'block';
        document.getElementById('analysis-descriptions').style.display = 'none';
        filterResults();
        renderPlotTabs();
        renderPlots();
    }

    window.filterResults = function() {
        selectedPhenotype = document.getElementById('phenotype-select').value;
        selectedGroupTest = document.getElementById('grouptest-select').value;
        selectedMaf = document.getElementById('maf-select').value;

        filteredData = allData.filter(row => 
            row.cohort === selectedCohort &&
            (selectedPhenotype === '' || row.phenotype === selectedPhenotype) &&
            (selectedGroupTest === '' || row.annot === selectedGroupTest) &&
            (selectedMaf === '' || row.max_maf === selectedMaf)
        );

        renderHitsTable(filteredData);
        renderPlotTabs();
        renderPlots();
    }

    function showDefaultView() {
        document.getElementById('default-view').style.display = 'block';
        document.getElementById('statistics').style.display = 'none';
        document.querySelector('.table-container').style.display = 'none';
        document.getElementById('plot-tabs').style.display = 'none';
        document.getElementById('cohort-buttons').style.display = 'none';
        document.getElementById('hits-table-title').style.display = 'none';
        document.getElementById('phenotype-filter').style.display = 'none';
        document.getElementById('analysis-descriptions').style.display = 'none';
    }

    function renderPlotTabs() {
        const plotTabs = document.getElementById('plot-tabs');
        const plotContainer = document.getElementById('plot-container');
        let activeTab = 'manhattan';

        function renderPlot() {
            plotTabs.innerHTML = `
                <div class="button-container">
                    <button id="manhattan-btn" class="button ${activeTab === 'manhattan' ? 'active' : ''}">Manhattan Plot</button>
                    <button id="qq-btn" class="button ${activeTab === 'qq' ? 'active' : ''}">QQ Plot</button>
                </div>
            `;
            plotContainer.innerHTML = `
                <img src="${activeTab === 'manhattan' ? `./Plots/${selectedCohort}.${selectedPhenotype}.manhattan_vertical.png` : `./Plots/${selectedCohort}.${selectedPhenotype}.qq.png`}" alt="${activeTab === 'manhattan' ? 'Manhattan Plot' : 'QQ Plot'}" class="mx-auto plot-image ${activeTab === 'manhattan' ? 'rotated-img' : ''}"/>
            `;

            document.getElementById('manhattan-btn').addEventListener('click', () => {
                activeTab = 'manhattan';
                renderPlot();
            });

            document.getElementById('qq-btn').addEventListener('click', () => {
                activeTab = 'qq';
                renderPlot();
            });
        }

        renderPlot();
    }

    function parseCSV(csv, delimiter = ',') {
        const lines = csv.split('\n');
        const headers = lines[0].split(delimiter).map(header => header.trim());
        const data = lines.slice(1).map(line => {
            const values = line.split(delimiter);
            return headers.reduce((obj, header, index) => {
                obj[header] = values[index] ? values[index].trim() : '';
                return obj;
            }, {});
        });
        return { headers, data };
    }
    

    function renderTableBody(data) {
        const tbody = document.querySelector('#hits-table tbody');
        const rowsPerPage = 10;
        const paginatedData = data.slice(0, rowsPerPage);
        
        tbody.innerHTML = paginatedData.map(row => `
            <tr>
                ${columns.map(col => `<td>${row[col] || ''}</td>`).join('')}
            </tr>
        `).join('');

        updatePagination(data.length);
    }

    function renderHitsTable(data) {
        const hitsTable = document.getElementById('hits-table');
        const headerRow = document.getElementById('header-row');
        const filterRow = document.getElementById('filter-row');

        headerRow.innerHTML = '';
        filterRow.innerHTML = '';

        columns.forEach(col => {
            headerRow.innerHTML += `<th>${col}</th>`;
            const input = document.createElement('input');
            input.type = 'text';
            input.placeholder = `Search ${col}`;
            input.dataset.column = col;
            input.value = currentFilters[col] || '';
            
            input.addEventListener('input', function() {
                currentFilters[col] = this.value;
                filterTable();
            });

            const th = document.createElement('th');
            th.appendChild(input);
            filterRow.appendChild(th);
        });

        renderTableBody(data);

        let currentPage = 0;
        const rowsPerPage = 10;

        function paginate(dataToRender) {
            const paginatedData = dataToRender.slice(currentPage * rowsPerPage, (currentPage + 1) * rowsPerPage);
            renderTableBody(paginatedData);

            document.getElementById('page-info').innerText = `Page ${currentPage + 1} of ${Math.ceil(dataToRender.length / rowsPerPage)}`;

            document.getElementById('first-page').disabled = currentPage === 0;
            document.getElementById('prev-page').disabled = currentPage === 0;
            document.getElementById('next-page').disabled = currentPage === Math.ceil(dataToRender.length / rowsPerPage) - 1;
            document.getElementById('last-page').disabled = currentPage === Math.ceil(dataToRender.length / rowsPerPage) - 1;
        }

        document.getElementById('first-page').addEventListener('click', () => {
            currentPage = 0;
            paginate(filteredData);
        });

        document.getElementById('prev-page').addEventListener('click', () => {
            if (currentPage > 0) {
                currentPage--;
                paginate(filteredData);
            }
        });

        document.getElementById('next-page').addEventListener('click', () => {
            if (currentPage < Math.ceil(filteredData.length / rowsPerPage) - 1) {
                currentPage++;
                paginate(filteredData);
            }
        });

        document.getElementById('last-page').addEventListener('click', () => {
            currentPage = Math.ceil(filteredData.length / rowsPerPage) - 1;
            paginate(filteredData);
        });

        paginate(filteredData);
    }

    function loadCSVFile() {
        console.log('Starting to load CSV file...');
        fetch('http://localhost:8000/saige_exwas_suggestive_regions.csv')
            .then(response => {
                console.log('CSV file fetched successfully', response);
                return response.text();
            })
            .then(csvText => {
                console.log('Parsing CSV text...');
                const { headers, data } = parseCSV(csvText);
                console.log('CSV parsed. Headers:', headers);
                columns = headers;
                allData = data;
                cohorts = [...new Set(allData.map(row => row.cohort))];
                console.log('Unique Cohorts:', cohorts);
                renderSideMenu();
                updateFilterOptions();
                showDefaultView(); // Show default view initially
            })
            .catch(error => {
                console.error('Error fetching the CSV file:', error);
                // Display an error message to the user
                document.getElementById('results-filter-submenu').innerHTML = '<p>Error loading cohorts. Please try again later.</p>';
            });
    }
    
    

    function renderSideMenu() {
        const submenu = document.getElementById('results-filter-submenu');
        submenu.innerHTML = '<h3>Cohorts:</h3>';
        cohorts.forEach(cohort => {
            const cohortLink = document.createElement('a');
            cohortLink.href = '#';
            cohortLink.textContent = cohort;
            cohortLink.onclick = (e) => {
                e.preventDefault();
                showCohort(cohort);
            };
            submenu.appendChild(cohortLink);
        });
        // Make the submenu visible
        submenu.style.display = 'block';
    }
    

    function updatePagination(totalRows) {
        const rowsPerPage = 10;
        const totalPages = Math.ceil(totalRows / rowsPerPage);
        const currentPage = 1;

        document.getElementById('page-info').innerText = `Page ${currentPage} of ${totalPages}`;

        document.getElementById('first-page').disabled = currentPage === 1;
        document.getElementById('prev-page').disabled = currentPage === 1;
        document.getElementById('next-page').disabled = currentPage === totalPages;
        document.getElementById('last-page').disabled = currentPage === totalPages;
    }

    function loadAnalysisDescriptions() {
        fetch('analysis_descriptions.json')
            .then(response => response.json())
            .then(data => {
                const container = document.getElementById('cards-container');
                container.innerHTML = '';
                data.forEach(description => {
                    const card = document.createElement('div');
                    card.className = 'card';
                    card.innerHTML = `
                        <h3>Cohorts: ${description.cohort}</h3>
                        <p>Mean: ${description.mean}</p>
                        <p>STD: ${description.std}</p>
                        <p>Min: ${description.min}</p>
                        <p>25th %ile: ${description.percentile_25}</p>
                        <p>Median: ${description.median}</p>
                        <p>75th %ile: ${description.percentile_75}</p>
                        <p>Max: ${description.max}</p>
                    `;
                    container.appendChild(card);
                });
            })
            .catch(error => console.error('Error fetching the JSON file:', error));
    }

    window.filterCards = function() {
        const filterValue = document.getElementById('filter-input').value.toUpperCase();
        const cards = document.querySelectorAll('#cards-container .card');
        cards.forEach(card => {
            const cohort = card.querySelector('h3').textContent.toUpperCase();
            if (cohort.includes(filterValue)) {
                card.style.display = 'block';
            } else {
                card.style.display = 'none';
            }
        });
    }

    function loadPlotManifest() {
        fetch('plots_manifest.csv')
            .then(response => response.text())
            .then(csvText => {
                const { headers, data } = parseCSV(csvText);
                plotManifest = data;
                updateFilterOptions();
            })
            .catch(error => console.error('Error fetching the plot manifest CSV file:', error));
    }

    function updateFilterOptions() {
        const phenotypes = [...new Set(allData.map(row => row.phenotype))];
        const groupTests = [...new Set(allData.map(row => row.annot))];
        const mafs = [...new Set(allData.map(row => row.max_maf))];
    
        renderFilterDropdown('phenotype-select', phenotypes);
        renderFilterDropdown('grouptest-select', groupTests);
        renderFilterDropdown('maf-select', mafs);
    }
    

    function renderFilterDropdown(id, options) {
        const select = document.getElementById(id);
        select.innerHTML = `<option value="">All</option>`;
        options.forEach(option => {
            select.innerHTML += `<option value="${option}">${option}</option>`;
        });
    }

    function renderPlots() {
        const plotContainer = document.getElementById('plot-container');
        plotContainer.innerHTML = '';

        const relevantPlots = plotManifest.filter(plot => 
            plot.cohort === selectedCohort &&
            plot.phenotype === selectedPhenotype &&
            plot.grouptest === selectedGroupTest &&
            plot.maf === selectedMaf &&
            plot.plots_file_present === 'True'
        );

        relevantPlots.forEach(plot => {
            const img = document.createElement('img');
            img.src = plot.expected_plots_file;
            img.alt = `${plot.plot_types} plot`;
            img.className = 'plot-image';
            plotContainer.appendChild(img);
        });

        if (relevantPlots.length === 0) {
            plotContainer.innerHTML = '<p>No plots available for the selected filters.</p>';
        }
    }

    // showDefaultView();
    loadCSVFile();
    loadAnalysisDescriptions();
    loadPlotManifest();
});