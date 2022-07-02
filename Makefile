test:
	@echo "    Running unit-tests"
	pytest tests

lint:
	@echo "    Linting tomplot codebase"
	@python3 -m flake8 tomplot
	@echo "    Linting tomplot tests"
	@python3 -m flake8 tests
