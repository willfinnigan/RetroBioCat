from retrobiocat_web.app.app import create_app
import os

production_mode = os.environ.get('PRODUCTION') or False

main_app = create_app(use_talisman=production_mode)

if __name__ == '__main__':
    main_app.run(debug=not production_mode)
