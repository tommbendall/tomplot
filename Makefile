test:
	@echo "    Running unit-tests"
	pytest tests

show-test:
	@echo "    Running unit-tests"
	pytest tests --show_plots

save-test:
	@echo "    Running unit-tests"
	pytest tests --save_plots

overwrite-test:
	@echo "    Running unit-tests"
	pytest tests --overwrite_plots

lint:
	@echo "    Linting tomplot tests"
	@python3 -m flake8 tests
	@echo "    Linting tomplot codebase"
	@python3 -m flake8 tomplot
